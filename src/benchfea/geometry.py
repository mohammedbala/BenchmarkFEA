"""Member local axes and flexible-length geometry (unit-agnostic).

Rule (docs/CONVENTIONS.md, section "Member local axes"):

* local x = unit vector from the flexible start (end i) to the flexible end (end j);
* reference vector v = ``orientation.referenceVector`` if given, otherwise the
  global vertical axis, except that for a vertical member (|x cross vertical|
  <= 1e-6) v = global +X;
* local z = unit(x cross v), local y = z cross x  (so y lies in the plane of x
  and v, on v's side);
* finally y and z are rotated about +x by the roll angle (right-hand rule).

The triad is returned as a 3x3 matrix R whose ROWS are the local x, y, z unit
vectors in global components:  v_local = R @ v_global.
"""
from __future__ import annotations

import math

import numpy as np

from .defaults import VERTICAL_DETECTION_SINE

AXIS = {"X": np.array([1.0, 0.0, 0.0]), "Y": np.array([0.0, 1.0, 0.0]), "Z": np.array([0.0, 0.0, 1.0])}


class GeometryError(ValueError):
    """Raised with a ``code`` telling validation what is wrong."""

    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def is_vertical(x_unit: np.ndarray, vertical: str) -> bool:
    return float(np.linalg.norm(np.cross(x_unit, AXIS[vertical]))) <= VERTICAL_DETECTION_SINE


def local_triad(p_i, p_j, vertical: str = "Z", reference_vector=None, roll_deg: float = 0.0,
                zero_length: float = 0.0) -> np.ndarray:
    """Local axes of a member running from ``p_i`` to ``p_j`` (see module doc)."""
    d = np.asarray(p_j, float) - np.asarray(p_i, float)
    length = float(np.linalg.norm(d))
    if not length > zero_length:
        raise GeometryError("member.zero-length", "member has zero length")
    x = d / length
    if reference_vector is not None:
        v = np.asarray(reference_vector, float)
        nv = float(np.linalg.norm(v))
        if nv == 0.0:
            raise GeometryError("member.orientation-collinear", "orientation.referenceVector is the zero vector")
        v = v / nv
        if float(np.linalg.norm(np.cross(x, v))) <= VERTICAL_DETECTION_SINE:
            raise GeometryError("member.orientation-collinear",
                                "orientation.referenceVector is parallel to the member axis")
    elif is_vertical(x, vertical):
        v = AXIS["X"]
    else:
        v = AXIS[vertical]
    z = np.cross(x, v)
    z /= np.linalg.norm(z)
    y = np.cross(z, x)
    if roll_deg:
        c, s = math.cos(math.radians(roll_deg)), math.sin(math.radians(roll_deg))
        y, z = c * y + s * z, -s * y + c * z
    return np.vstack([x, y, z])


def member_geometry(member: dict, nodes: dict, vertical: str, zero_length: float = 0.0):
    """Geometry of one member from raw document data (any consistent length unit).

    Returns ``(start, end, R, length, offset_i, offset_j)``: the flexible-length
    end points, the local triad (rows), the flexible length, and the rigid end
    offsets as global vectors.  Offsets given in local axes are interpreted in
    the triad of the node-to-node line.  Raises :class:`GeometryError`.
    """
    def xyz(n):
        return np.array([n["x"], n["y"], n.get("z", 0.0)], float)

    p_i, p_j = xyz(nodes[member["i"]]), xyz(nodes[member["j"]])
    ori = member.get("orientation") or {}
    ref, roll = ori.get("referenceVector"), ori.get("roll", 0.0)
    off = member.get("offsets") or {}
    o_i = np.asarray(off.get("i", (0.0, 0.0, 0.0)), float)
    o_j = np.asarray(off.get("j", (0.0, 0.0, 0.0)), float)
    if off.get("axes", "global") == "local" and (o_i.any() or o_j.any()):
        try:
            Rn = local_triad(p_i, p_j, vertical, ref, roll, zero_length)
        except GeometryError as e:
            if e.code == "member.zero-length":
                raise GeometryError("member.zero-length", "offsets in local axes need distinct end nodes") from None
            raise
        o_i, o_j = Rn.T @ o_i, Rn.T @ o_j
    start, end = p_i + o_i, p_j + o_j
    try:
        R = local_triad(start, end, vertical, ref, roll, zero_length)
    except GeometryError as e:
        if e.code == "member.zero-length" and (o_i.any() or o_j.any()):
            raise GeometryError("member.offset-zero-length",
                                "the rigid end offsets leave no flexible length") from None
        raise
    return start, end, R, float(np.linalg.norm(end - start)), o_i, o_j


def model_size(nodes: dict) -> float:
    """Diagonal of the bounding box of all nodes (0 if not computable)."""
    pts = np.array([[n["x"], n["y"], n.get("z", 0.0)] for n in nodes.values()], float)
    return float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0))) if len(pts) else 0.0
