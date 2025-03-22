#! /usr/bin/env python3
#
# KiCad part/footprint importer
#
import os
import glob

from math import *
import shapely.geometry as sg
import shapely.affinity as sa
import shapely.ops as so

from .sexp_parser import *
from pcbflow import *

KI_LAYER_DICT = {
    "F.SilkS": "GTO",
    "F.Paste": "GTP",
    "F.Mask": "GTS",
    "F.Cu": "GTL",
    "F.Fab": "GTD",
}

FP_LIB_PATH = None
ALL_KICAD_MOD_FILES = None


class KiCadPart(PCBPart):
    def __init__(self, dc, val=None, source=None, libraryfile=None, **kwargs):
        self.libraryfile = libraryfile
        self.smd_pads = []
        self.pin_pads = []
        self.polys = []
        self.lines = []
        self.circles = []
        self.labels = []
        self.docu = []
        if "side" not in self.__dict__:
            self.side = "top"
        if "side" in kwargs:
            self.side = kwargs["side"]
        self.parse()
        if "family" in kwargs:
            self.family = kwargs["family"]
        super().__init__(dc, val, source, **kwargs)

    def _add_obj_to_layer(self, obj, layer):
        if layer == "GTO":
            self.board.get_silk_layer(side=self.side).add(obj)
        elif layer == "GTD":
            self.board.get_docu_layer(side=self.side).add(obj)
        elif layer == "GTP":
            self.board.get_paste_layer(side=self.side).add(obj)
        elif layer == "GTL":
            if self.side == "bottom":
                self.board.layers["GBL"].add(obj)
            else:
                self.board.layers["GTL"].add(obj)

    def place(self, dc):
        for line in self.lines:
            p0 = dc.copy().goxy(*line["coords"][0])
            p1 = dc.copy().goxy(*line["coords"][1])
            width = self.board.drc.silk_width
            if line["width"] > 0:
                width = line["width"]
            if width==0:
                print("warning: line width is zero")
            g = sg.LineString([p0.xy, p1.xy]).buffer(width / 2)
            self._add_obj_to_layer(g, line["layers"][0])

        for poly in self.polys:
            width = self.board.drc.silk_width
            if poly["width"] > 0:
                width = poly["width"]
            coords = []
            xyc = self.center.xy
            for c in poly["coords"]:
                coords.append((xyc[0] + c[0], xyc[1] + c[1]))
            if poly['fill']:
                g = sg.Polygon(coords).buffer(width/2)
            else:
                g = sg.Polygon(coords).exterior.buffer(width/2)
            self._add_obj_to_layer(g, poly["layer"])

        for circle in self.circles:
            width = self.board.drc.silk_width
            if circle["width"] > 0:
                width = circle["width"]
            xyc=self.center.copy().goxy(*circle["center"]).xy
            gc = sg.Point(xyc).buffer(circle["diameter"] / 2)
            if circle['fill']:
                g = sg.Polygon(gc.exterior.coords).buffer(width/2)
            else:
                g = gc.exterior.buffer(width/2)
            self._add_obj_to_layer(g, circle["layer"])

        for pad in self.smd_pads:
            p = dc.copy().goxy(*pad["xy"])
            p.rect(*pad["size"])
            p.set_name(pad["name"])
            if "GTL" in pad["layers"]:
                no_paste = True if "GTP" not in pad["layers"] else False
                self.smd_pad(p, ignore_paste=no_paste)
            elif "GTP" in pad["layers"]:
                self.board.get_paste_layer(side=self.side).add(p.poly())

        for pad in self.pin_pads:
            diameter = pad["size"][0]
            dc.push()
            dc.goxy(*pad["xy"])
            dc.board.add_drill(dc.xy, pad["drill"])
            shape = pad["shape"]
            if shape in ["long", "circle", "octagon", "rect"]:
                n = {"long": 60, "circle": 60, "octagon": 8, "rect": 4}[shape]
                if shape == "rect":
                    diameter /= 1.1
                p = dc.copy()
                p.n_agon(diameter / 2, n)
                p.set_name(pad["name"])
                p.part = self.id
                self.pads.append(p)
                p.pin_pad()
                dc.pop()

        if len(self.labels) > 0:
            for label in self.labels:
                #xyc = self.center.xy
                xy=self.center.copy().goxy(*label["xy"]).xy
                #xy = (xyc[0] + label["xy"][0], xyc[1] + label["xy"][1])
                self.board.add_text(
                    xy,
                    self.id,
                    angle=0,
                    scale=1.0,
                    side=self.side,
                    justify="center",
                )

    def _map_layers(self, layers):
        # Ensure 'layers' is always a list
        if isinstance(layers, str):
            layers = [layers]  # Convert single string into a list

        ml = []
        for layer in layers:
            clean_layer = layer.strip('"')  # Remove any extra quotes
            if clean_layer in KI_LAYER_DICT:
                ml.append(KI_LAYER_DICT[clean_layer])

        return ml

    def _parse_fp_text(self, items):
        xy = []
        layer = "F.Cu"  # Default layer if mapping fails
        text = items[0]

        for e in items:
            if isinstance(e, dict):
                if "at" in e:
                    xy = float(e["at"][0]), -float(e["at"][1])
                elif "layer" in e:
                    mapped_layers = self._map_layers(e["layer"])
                    if mapped_layers:
                        layer = mapped_layers[0]  # Use the first mapped layer
                    else:
                        print(f"Warning: Could not map layer '{e['layer']}', using default '{layer}'")

        if text == "reference":
            self.labels.append({"xy": xy, "text": text, "layer": layer})

    def _parse_fp_poly(self, items):
        coords = []
        fill = True
        width = 0.0  # Default value
        layer = "F.Cu"  # Default layer (or use a sensible default)

        for e in items:
            if isinstance(e, dict):
                if "pts" in e:
                    for pt in e["pts"]:
                        coords.append((float(pt[2]), -float(pt[3])))

                elif "width" in e:
                    width = float(e["width"][0])

                elif "stroke" in e:
                    for ln, key, val in e["stroke"]:
                        if "width" in key:
                            width = float(val)

                elif "fill" in e:
                    if "none" in e["fill"]:
                        fill = False

                elif "layer" in e:
                    mapped_layers = self._map_layers(e["layer"])
                    if mapped_layers:  # Ensure the list is not empty
                        layer = mapped_layers[0]
                    else:
                        print(f"Warning: Unknown layer '{e['layer']}', using default '{layer}'")

        self.polys.append({"coords": coords, "width": width, "layer": layer, "fill": fill})

    def _parse_fp_circle(self, items):
        center = (0, 0)
        width = 0
        diameter = 0
        layers = []
        fill=True
        for e in items:
            if isinstance(e, dict):
                if "center" in e:
                    pt = e["center"]
                    center = (float(pt[0]), -float(pt[1]))
                if "end" in e:
                    pt = e["end"]
                    end = (float(pt[0]), -float(pt[1]))
                    diameter = 2*hypot(center[0] - end[0], center[1] - end[1])
                elif "width" in e:
                    width = float(e["width"][0])
                elif "stroke" in e:
                    for ln,key,val in e['stroke']:
                        if key=="width":
                            width=float(val)
                elif "fill" in e:
                    if 'none' in e['fill']:
                        fill=False
                elif "layer" in e:
                    layer = self._map_layers(e["layer"])[0]
        self.circles.append(
            {"center": center, "diameter": diameter, "width": width, "layer": layer, "fill" : fill}
        )

    def _parse_fp_line(self, items):
        coords = []
        width = 0
        layers = []
        for e in items:
            if isinstance(e, dict):
                if "start" in e:
                    coord = e["start"]
                    x, y = float(coord[0]), float(coord[1])
                    coords.append((float(coord[0]), -float(coord[1])))
                elif "end" in e:
                    coord = e["end"]
                    x, y = float(coord[0]), float(coord[1])
                    coords.append((float(coord[0]), -float(coord[1])))
                elif "width" in e:
                    width = float(e["width"][0])
                elif "stroke" in e:
                    for ln,key,val in e['stroke']:
                        if key=="width":
                            width=float(val)
                elif "layer" in e:
                    layers = self._map_layers(e["layer"])
        if len(layers) > 0:
            self.lines.append({"coords": coords, "width": width, "layers": layers})

    def _parse_pad(self, items):
        xy = (0, 0)
        size = (0, 0)
        name = items[0]
        layers = []
        shape = ""
        drill = 0
        if items[1] == "smd":
            for e in items[2:]:
                if isinstance(e, dict):
                    if "at" in e:
                        xy = float(e["at"][0]), -float(e["at"][1])
                    elif "size" in e:
                        size = float(e["size"][0]), float(e["size"][1])
                    elif "layers" in e:
                        layers = self._map_layers(e["layers"])
            self.smd_pads.append(
                {"name": name, "xy": xy, "size": size, "layers": layers}
            )

        elif items[1] == "thru_hole":
            shape = items[2]
            for e in items[3:]:
                if isinstance(e, dict):
                    if "at" in e:
                        xy = float(e["at"][0]), -float(e["at"][1])
                    elif "size" in e:
                        size = float(e["size"][0]), float(e["size"][1])
                    # TODO: support oval shaped drill, i.e. slot
                    elif "drill" in e:
                        try:
                            drill = float(e["drill"][0])
                        except:
                            pass
            self.pin_pads.append(
                {"name": name, "xy": xy, "size": size, "drill": drill, "shape": shape}
            )
        elif items[1] == "np_thru_hole":
            shape = items[2]
            for e in items[3:]:
                if isinstance(e, dict):
                    if "at" in e:
                        xy = float(e["at"][0]), -float(e["at"][1])
                    elif "size" in e:
                        size = float(e["size"][0]), float(e["size"][1])
                    # TODO: support oval shaped drill, i.e. slot
                    elif "drill" in e:
                        try:
                            drill = float(e["drill"][0])
                        except:
                            pass
            self.pin_pads.append(
                {"name": name, "xy": xy, "size": size, "drill": drill, "shape": shape}
            )

    def parse(self):
        with open(self.libraryfile, "r") as f:
            module = parseSexp(f.read())
        for idx, item in enumerate(module):
            if not isinstance(item, list):
                if item == "module":
                    self.footprint = module[idx + 1]
                    self.family = infer_family(self.footprint)
                    break

        d = {str(k[0]): k[1:] for k in module if isinstance(k, list)}
        md = {}
        for k, v in d.items():
            vd = []
            for kk in v[1:]:
                if isinstance(kk, list):
                    vd.append({kk[1]: kk[2:]})
                else:
                    vd.append(kk)
            md[k] = {v[0]: vd}

        for k, v in md.items():
            if "fp_line" in v:
                self._parse_fp_line(v["fp_line"])
            if "pad" in v:
                self._parse_pad(v["pad"])
            if "fp_text" in v:
                self._parse_fp_text(v["fp_text"])
            if "fp_poly" in v:
                self._parse_fp_poly(v["fp_poly"])
            if "fp_circle" in v:
                self._parse_fp_circle(v["fp_circle"])


# TODO   (fp_arc (start 0 0) (end 0 4) (angle -65) (layer F.Fab) (width 0.1))


class SkiPart(KiCadPart):
    def __init__(self, dc, skipart, **kwargs):
        self.skipart = skipart
        self.kicad_path = os.environ["KICAD_DIR"]
        if not self.kicad_path:
            raise FileNotFoundError("KiCad directory not found. Please set environment variable: KICAD_DIR")
        self.footprint = skipart.footprint
        self.libraryfile = self._find_footprint_file(skipart.footprint)
        super().__init__(
            dc,
            val=skipart.value,
            libraryfile=self.libraryfile,
            family=skipart.ref_prefix,
            ref=skipart.ref,
            **kwargs
        )
        self._assign_pads()

    def _find_footprint_file(self, footprint_name):
        """
        Find the footprint file associated with the given footprint name.
        Modern KiCad versions (v5+) use library files and paths in a more structured way.
        """

        footprint_path = self._search_footprint_in_project(footprint_name)
        if footprint_path:
            return footprint_path

        raise FileNotFoundError(f"Footprint file for '{footprint_name}' not found.")

    def _search_footprint_in_project(self, footprint):
        """
        Search for the footprint within the current project directory.
        """

        kicad_dir = self.kicad_path
        footprint_dir = os.path.join(kicad_dir, "footprints")

        footprint_files = glob.glob(os.path.join(footprint_dir, f"**/{footprint}*.kicad_mod"), recursive=True)

        if footprint_files:
            return footprint_files[0] # TODO resolve if multiple footprints with same names exist
        else:
            return None

    def _assign_pads(self):
        """
        Assign the correct nets to the pads based on the schematic pin connections.
        This part of the logic remains the same, where pads are matched to their respective pins.
        """
        for pad in self.pads:
            for pin in self.skipart.pins:
                if str(pin.num) == str(pad.name):
                    if len(pin.nets) == 1:
                        pad.name = pin.nets[0].name
