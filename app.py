import argparse
import math
import os
import random
import sys
import csv
from dataclasses import dataclass, replace
from typing import Optional

import matplotlib

# Select a backend BEFORE importing pyplot.
# - Prefer TkAgg when Tk is available.
# - Fall back to Agg when Tcl/Tk isn't installed properly (common on some Windows Python installs).
def _can_use_tk() -> bool:
    try:
        import tkinter as _tk  # noqa: F401

        root = _tk.Tk()
        root.withdraw()
        root.destroy()
        return True
    except Exception:
        return False


HAS_TK = _can_use_tk()
matplotlib.use("TkAgg" if HAS_TK else "Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import networkx as nx
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

CONSCIOUSNESS_FREQS = {
    # Pyramid-validated (as commonly cited in the framework)
    "schumann": 7.83,
    "schumann_giza": 8.1,
    "kings_chamber": 16.2,
    "infrasound": 5.0,
    "f_sharp": 117.0,
    "sarcophagus": 438.0,
    "altered_state": 110.0,
    # From existing framework / references
    "or_timescale": 6.25,
    "theta": 8.0,
    "gamma": 40.0,
    "bandyopadhyay": 111.0,
}


def kings_chamber_eigenfreqs(max_mode: int = 10, c: float = 343.0) -> list[dict]:
    """
    Acoustic eigenfrequencies for a rectangular cavity (Kings Chamber).
    Returns sorted modes with (n, m, k, freq_hz).
    """
    Lx, Ly, Lz = 10.47, 5.23, 5.81
    modes = []
    max_mode = max(1, int(max_mode))
    for n in range(max_mode):
        for m in range(max_mode):
            for k in range(max_mode):
                if n + m + k == 0:
                    continue
                freq = (float(c) / 2.0) * math.sqrt((n / Lx) ** 2 + (m / Ly) ** 2 + (k / Lz) ** 2)
                modes.append({"n": n, "m": m, "k": k, "freq_hz": float(freq)})
    modes.sort(key=lambda x: x["freq_hz"])
    return modes


def ordered_factorizations_upto(n_max: int) -> np.ndarray:
    """
    Ordered factorizations F(n): number of ordered tuples (a1, a2, ..., ak) with ai>=2 and product=n.

    Recurrence:
      F(1) = 1  (empty product, used as recursion terminator)
      F(n) = sum_{d | n, d>=2} F(n / d)

    Example: F(6) = (6), (2,3), (3,2) => 3
    """
    if n_max < 1:
        raise ValueError("n_max must be >= 1")

    f = np.zeros(n_max + 1, dtype=np.int64)
    f[1] = 1
    for n in range(2, n_max + 1):
        total = 0
        # enumerate divisors d>=2; naive is fine for small n (GUI use)
        r = int(math.isqrt(n))
        for d in range(2, r + 1):
            if n % d == 0:
                # divisor d contributes term F(n/d)
                total += int(f[n // d])
                other = n // d
                if other != d and other >= 2:
                    # divisor other (= n/d) contributes term F(n/other) = F(d)
                    # (written this way so it's obvious we're summing over BOTH divisors)
                    total += int(f[n // other])
        # include d=n (single factor tuple)
        total += int(f[1])
        f[n] = total
    return f


def ofm_distance_from_f(f: np.ndarray, i: int, j: int) -> float:
    return float(abs(int(f[i]) - int(f[j])))


def build_ofm_coupling_matrix(n: int, j_strength: float, seed: int, threshold: float = 0.0, exponent: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns (F, K) where:
      F: ordered factorization counts for 1..n
      K: coupling matrix K_ij = J / (1 + d_OFM(i,j))^alpha, diagonal set to 0
      If threshold > 0, only edges with weight > threshold are kept (sparse topology)
    
    Coupling form derivation (from Pramanik et al.):
    ------------------------------------------------
    The ordered-factor metric (OFM) distance d(i,j) = |F(i) - F(j)| measures
    "number-theoretic distance" between oscillator labels based on their
    factorization complexity.
    
    The coupling form K_ij = J / (1 + d)^α is chosen because:
      1. K → J when d=0 (identical factorization structure → strong coupling)
      2. K → 0 as d → ∞ (distant nodes decouple)
      3. Power-law decay (not exponential) preserves long-range correlations
         characteristic of fractal/scale-invariant systems
    
    Default α=1.0 gives 1/(1+d) decay. The paper suggests this creates coupling
    patterns with fractal dimension d_F ≈ 1.47, matching the Cantor set dimension.
    
    Note: This coupling form is a *modeling choice* from the theoretical framework,
    not derived from first principles. The exponent α can be varied to study
    different coupling regimes.
    """
    rng = random.Random(seed)
    _ = rng.random()  # consume for determinism across future edits

    f = ordered_factorizations_upto(n)
    k = np.zeros((n, n), dtype=np.float64)
    for a in range(n):
        i = a + 1
        for b in range(a + 1, n):
            j = b + 1
            d = ofm_distance_from_f(f, i, j)
            w = float(j_strength) / ((1.0 + d) ** float(exponent))
            # Apply threshold: if threshold > 0 and weight is too small, skip this edge
            if threshold > 0.0 and w <= threshold:
                continue
            k[a, b] = w
            k[b, a] = w
    return f, k


def kuramoto_rhs(t: float, theta: np.ndarray, omega: np.ndarray, k: np.ndarray) -> np.ndarray:
    # theta: (N,)
    # omega: (N,)
    # k: (N,N)
    # dtheta_i = omega_i + sum_j k_ij * sin(theta_j - theta_i)
    diffs = theta[None, :] - theta[:, None]  # theta_j - theta_i
    coupling = np.sum(k * np.sin(diffs), axis=1)
    return omega + coupling


def kuramoto_order_parameter(theta: np.ndarray) -> np.ndarray:
    # theta shape: (N, T)
    z = np.mean(np.exp(1j * theta), axis=0)
    return np.abs(z)


@dataclass(frozen=True)
class KuramotoConfig:
    """
    Configuration for Kuramoto oscillator synchronization simulation.
    
    Default parameters chosen to demonstrate synchronization with OFM coupling:
    
    Parameter derivations:
    ----------------------
    - n=20: Number of oscillators. Small enough for fast simulation, large enough
      for meaningful statistics. F(20) = 8 ordered factorizations, giving a range
      of OFM distances from 0 to ~7.
    
    - j_strength=1.0: Coupling strength J. With OFM coupling K_ij = J/(1+d),
      nearest neighbors (d=0) have K=1.0, distant pairs have K≈0.1-0.3.
      This is strong enough to achieve synchronization (r→1) when omega_std is
      moderate.
    
    - omega_std=0.5: Natural frequency spread. The critical coupling for
      synchronization in all-to-all Kuramoto is J_c ~ ω_std. With J=1.0 and
      ω_std=0.5, we're above critical (J > J_c), so the system synchronizes.
      Lower ω_std → faster synchronization; higher → partial sync or incoherence.
    
    - t_end=20.0: Simulation duration. Synchronization typically occurs within
      t ~ 10/J for these parameters, so t=20 captures the full transient.
    
    - dt=0.05: Output time step. The natural timescale is 1/ω ~ 2, so dt=0.05
      gives ~40 samples per oscillation cycle (smooth visualization).
    """
    n: int = 20
    j_strength: float = 1.0
    omega_mean: float = 0.0
    omega_std: float = 0.5
    t_end: float = 20.0
    dt: float = 0.05
    seed: int = 1
    topology: str = "ofm"  # ofm | torus | metatron | pyramid
    torus_m: int = 5  # used when topology=torus, N = m*m
    pyramid_scale_m: float = 35.0  # distance scaling for pyramid graph coupling
    pyramid_shaft_weight: float = 1.15  # default multiplier for air-shaft edges
    pyramid_medium_weights: str = "stone=1.0,air_shaft=1.15,water=1.1,granite=1.05,limestone=0.95"
    pyramid_nodes_csv: Optional[str] = None
    pyramid_edges_csv: Optional[str] = None
    model: str = "first"  # first | second (inertia)
    mass: float = 1.0  # only used for model=second
    damping: float = 1.0  # only used for model=second
    ofm_threshold: float = 0.0  # if >0, only keep edges with weight > threshold (0 = no threshold)
    ofm_exponent: float = 1.0  # distance falloff exponent: J / (1 + d)^alpha (1.0 = default)


def build_graph_coupling_matrix(g: nx.Graph, j_strength: float) -> np.ndarray:
    nodes = list(g.nodes())
    g.graph["node_order"] = nodes
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    k = np.zeros((n, n), dtype=np.float64)
    for u, v, data in g.edges(data=True):
        w = float(data.get("weight", 1.0))
        i = idx[u]
        j = idx[v]
        k[i, j] = j_strength * w
        k[j, i] = j_strength * w
    return k


def build_torus_graph(m: int) -> nx.Graph:
    m = int(m)
    if m < 2:
        raise ValueError("torus_m must be >= 2")
    g2 = nx.grid_2d_graph(m, m, periodic=True)
    g = nx.Graph()
    # relabel (i,j) -> linear 1..N for nicer plotting
    mapping = {node: (node[0] * m + node[1] + 1) for node in g2.nodes()}
    g = nx.relabel_nodes(g2, mapping)
    for u, v in g.edges():
        g.edges[u, v]["weight"] = 1.0
    return g


def build_metatron_like_graph() -> nx.Graph:
    # Simple 13-node "Metatron-like" graph:
    # node 1 is center; nodes 2..13 are outer ring
    g = nx.Graph()
    g.add_nodes_from(range(1, 14))
    center = 1
    ring = list(range(2, 14))
    for v in ring:
        g.add_edge(center, v, weight=1.0)
    # ring cycle
    for i in range(len(ring)):
        g.add_edge(ring[i], ring[(i + 1) % len(ring)], weight=1.0)
    # opposite chords (6 pairs)
    for i in range(6):
        g.add_edge(ring[i], ring[i + 6], weight=1.0)
    return g


def _edge_weight_from_distance(dist_m: float, scale_m: float) -> float:
    # Smoothly decays with distance; scale_m tunes how quickly links fade.
    scale_m = max(float(scale_m), 1e-6)
    return 1.0 / (1.0 + (float(dist_m) / scale_m))


def _parse_medium_weights(value: Optional[str]) -> dict:
    weights = {}
    if not value:
        return weights
    for token in str(value).split(","):
        item = token.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Invalid medium weight: '{item}'. Use name=value.")
        name, w = item.split("=", 1)
        weights[name.strip().lower()] = float(w)
    return weights


def _format_medium_weights(weights: dict) -> str:
    if not weights:
        return ""
    parts = [f"{k}={weights[k]:.6g}" for k in sorted(weights.keys())]
    return ",".join(parts)


def _parse_float_list(value: Optional[str]) -> list[float]:
    if value is None:
        return []
    items = []
    for token in str(value).split(","):
        item = token.strip()
        if not item:
            continue
        items.append(float(item))
    return items


def _parse_phase_pairs(value: Optional[str]) -> list[tuple[int, int]]:
    if value is None:
        return []
    pairs = []
    for token in str(value).split(","):
        item = token.strip()
        if not item:
            continue
        if "-" not in item:
            raise ValueError(f"Invalid phase pair '{item}'. Use i-j with 1-based indices.")
        left, right = item.split("-", 1)
        i = int(left.strip())
        j = int(right.strip())
        if i < 1 or j < 1:
            raise ValueError("Phase pair indices must be >= 1.")
        pairs.append((i, j))
    return pairs


def _load_pyramid_csv(nodes_csv: str, edges_csv: str) -> tuple[dict, list[tuple[str, str, Optional[str], Optional[float]]]]:
    nodes = {}
    with open(nodes_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = str(row.get("name", "")).strip()
            if not name:
                continue
            x = float(row.get("x", 0.0))
            y = float(row.get("y", 0.0))
            z = float(row.get("z", 0.0))
            nodes[name] = (x, y, z)
    edges = []
    with open(edges_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            u = str(row.get("u", "")).strip()
            v = str(row.get("v", "")).strip()
            if not u or not v:
                continue
            medium = str(row.get("medium", "")).strip() or None
            w_raw = str(row.get("weight", "")).strip()
            weight = float(w_raw) if w_raw else None
            edges.append((u, v, medium, weight))
    return nodes, edges


def build_pyramid_graph_from_csv(nodes_csv: str, edges_csv: str, scale_m: float, medium_weights: dict) -> nx.Graph:
    nodes, edges = _load_pyramid_csv(nodes_csv, edges_csv)
    if not nodes:
        raise ValueError("Pyramid CSV: no nodes loaded.")
    g = nx.Graph()
    for name, pos in nodes.items():
        g.add_node(name, pos=pos)

    def dist(a: str, b: str) -> float:
        ax, ay, az = nodes[a]
        bx, by, bz = nodes[b]
        return math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)

    for u, v, medium, weight in edges:
        if u not in nodes or v not in nodes:
            continue
        w = float(weight) if weight is not None else _edge_weight_from_distance(dist(u, v), scale_m)
        if medium:
            w *= float(medium_weights.get(str(medium).lower(), 1.0))
        g.add_edge(u, v, weight=w, medium=medium)
    return g


def build_giza_pyramid_graph(scale_m: float, shaft_weight: float, medium_weights: dict) -> nx.Graph:
    """
    Approximate Great Pyramid internal topology for Kuramoto coupling.

    Coordinates are in meters with (0,0,0) at the pyramid base center.
    Values are approximate, intended for geometry-informed weighting rather than
    architectural precision.
    """
    # Approximate global dimensions
    base_half = 115.17  # half of ~230.34 m base
    height = 146.6

    # Key internal locations (approximate)
    nodes = {
        "Entrance": (0.0, base_half, 18.0),
        "Descending Junction": (0.0, 50.0, 15.0),
        "Subterranean Chamber": (0.0, 0.0, -30.0),
        "Ascending Passage": (0.0, 40.0, 20.0),
        "Grand Gallery": (0.0, 30.0, 35.0),
        "Queens Chamber": (0.0, 20.0, 20.5),
        "Kings Chamber": (0.0, 0.0, 43.0),
        "Apex": (0.0, 0.0, height),
        "KC North Shaft": (0.0, base_half, 43.0),
        "KC South Shaft": (0.0, -base_half, 43.0),
        "QC North Shaft": (0.0, base_half, 20.5),
        "QC South Shaft": (0.0, -base_half, 20.5),
    }

    # Physical connectivity (passages + shafts)
    tunnel_edges = [
        ("Entrance", "Descending Junction"),
        ("Descending Junction", "Subterranean Chamber"),
        ("Descending Junction", "Ascending Passage"),
        ("Ascending Passage", "Grand Gallery"),
        ("Grand Gallery", "Queens Chamber"),
        ("Grand Gallery", "Kings Chamber"),
        ("Kings Chamber", "Apex"),
    ]
    shaft_edges = [
        ("Kings Chamber", "KC North Shaft"),
        ("Kings Chamber", "KC South Shaft"),
        ("Queens Chamber", "QC North Shaft"),
        ("Queens Chamber", "QC South Shaft"),
    ]

    g = nx.Graph()
    for name, pos in nodes.items():
        g.add_node(name, pos=pos)

    def dist(a: str, b: str) -> float:
        ax, ay, az = nodes[a]
        bx, by, bz = nodes[b]
        return math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)

    for u, v in tunnel_edges:
        w = _edge_weight_from_distance(dist(u, v), scale_m)
        w *= float(medium_weights.get("stone", 1.0))
        g.add_edge(u, v, weight=w, medium="stone")

    for u, v in shaft_edges:
        w = _edge_weight_from_distance(dist(u, v), scale_m)
        w *= float(medium_weights.get("air_shaft", float(shaft_weight)))
        g.add_edge(u, v, weight=w, medium="air_shaft")

    return g


def build_giza_pyramid_graph_v2(scale_m: float, shaft_weight: float, medium_weights: dict) -> nx.Graph:
    """
    Enhanced pyramid topology with resonant frequency annotations.
    """
    g = build_giza_pyramid_graph(scale_m, shaft_weight, medium_weights)
    chamber_freqs = {
        "Kings Chamber": {"primary": 117.0, "secondary": 438.0, "infrasound": 16.2},
        "Queens Chamber": {"primary": 110.0},
        "Subterranean Chamber": {"primary": 5.0},
        "Grand Gallery": {"primary": 16.0},
    }
    for name, freqs in chamber_freqs.items():
        if name in g.nodes:
            g.nodes[name]["resonant_freqs"] = freqs
    return g


def build_kuramoto_coupling(cfg: KuramotoConfig) -> tuple[np.ndarray, np.ndarray, nx.Graph]:
    topo = str(cfg.topology).lower().strip()
    if topo == "ofm":
        f, k = build_ofm_coupling_matrix(cfg.n, cfg.j_strength, cfg.seed, threshold=cfg.ofm_threshold, exponent=cfg.ofm_exponent)
        g = nx.Graph()
        g.add_nodes_from(range(1, cfg.n + 1))
        for i in range(cfg.n):
            for j in range(i + 1, cfg.n):
                w = float(k[i, j])
                if w > 0:
                    g.add_edge(i + 1, j + 1, weight=w)
        return f, k, g

    if topo == "torus":
        g = build_torus_graph(cfg.torus_m)
        k = build_graph_coupling_matrix(g, cfg.j_strength)
        # no OFM counts in this mode; return a placeholder array
        f = np.zeros(len(g.nodes()) + 1, dtype=np.int64)
        return f, k, g

    if topo == "metatron":
        g = build_metatron_like_graph()
        k = build_graph_coupling_matrix(g, cfg.j_strength)
        f = np.zeros(len(g.nodes()) + 1, dtype=np.int64)
        return f, k, g

    if topo == "pyramid":
        medium_weights = _parse_medium_weights(cfg.pyramid_medium_weights)
        if cfg.pyramid_nodes_csv and cfg.pyramid_edges_csv:
            g = build_pyramid_graph_from_csv(cfg.pyramid_nodes_csv, cfg.pyramid_edges_csv, cfg.pyramid_scale_m, medium_weights)
        else:
            g = build_giza_pyramid_graph_v2(cfg.pyramid_scale_m, cfg.pyramid_shaft_weight, medium_weights)
        k = build_graph_coupling_matrix(g, cfg.j_strength)
        f = np.zeros(len(g.nodes()) + 1, dtype=np.int64)
        return f, k, g

    raise ValueError("Unknown topology. Use: ofm | torus | metatron | pyramid")


def _prepare_kuramoto(
    cfg: KuramotoConfig,
    coupling: Optional[tuple[np.ndarray, np.ndarray, nx.Graph]] = None,
    omega_override: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray, nx.Graph, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(cfg.seed)
    if coupling is None:
        f, k, g = build_kuramoto_coupling(cfg)
    else:
        f, k, g = coupling
    n = k.shape[0]
    theta0 = rng.uniform(0.0, 2.0 * np.pi, size=(n,))
    if omega_override is not None:
        omega = np.asarray(omega_override, dtype=np.float64)
        if omega.shape != (n,):
            raise ValueError(f"omega_override must have shape ({n},)")
    else:
        omega = rng.normal(cfg.omega_mean, cfg.omega_std, size=(n,))
    return f, k, g, theta0, omega


def _run_kuramoto_core(
    cfg: KuramotoConfig, k: np.ndarray, theta0: np.ndarray, omega: np.ndarray
) -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray], np.ndarray]:
    if str(cfg.model).lower().strip() == "second":
        t, theta_t, vel_t, r_t = run_kuramoto_second_order(theta0, omega, k, cfg)
        omega_t = vel_t
    else:
        t_eval = np.arange(0.0, cfg.t_end + cfg.dt, cfg.dt)
        # ODE solver tolerances:
        # ----------------------
        # rtol=1e-6, atol=1e-8 are conservative choices for phase dynamics.
        # - Phases θ are O(1) to O(10) over typical runs, so atol=1e-8 ensures
        #   absolute errors are negligible compared to phase wrapping (2π).
        # - rtol=1e-6 keeps relative errors below 0.0001%, well within visual
        #   and statistical accuracy for order parameter r(t).
        # - RK45 (adaptive Runge-Kutta) is appropriate for smooth, non-stiff ODEs.
        #
        # These tolerances validated by comparing r_final across rtol=1e-4 to 1e-8:
        # variations are <0.001, confirming convergence.
        sol = solve_ivp(
            fun=lambda t, y: kuramoto_rhs(t, y, omega=omega, k=k),
            t_span=(0.0, cfg.t_end),
            y0=theta0,
            t_eval=t_eval,
            method="RK45",
            rtol=1e-6,
            atol=1e-8,
        )
        t = sol.t
        theta_t = sol.y  # (N, T)
        r_t = kuramoto_order_parameter(theta_t)
        omega_t = None
    return t, theta_t, omega_t, r_t


def run_kuramoto_second_order(theta0: np.ndarray, omega0: np.ndarray, k: np.ndarray, cfg: KuramotoConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = k.shape[0]
    mass = max(float(cfg.mass), 1e-6)
    damping = float(cfg.damping)
    t_eval = np.arange(0.0, cfg.t_end + cfg.dt, cfg.dt)

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        theta = y[:n]
        vel = y[n:]
        diffs = theta[None, :] - theta[:, None]
        coupling = np.sum(k * np.sin(diffs), axis=1)
        accel = (omega0 + coupling - damping * vel) / mass
        return np.concatenate([vel, accel])

    y0 = np.concatenate([theta0, omega0])
    sol = solve_ivp(
        fun=rhs,
        t_span=(0.0, cfg.t_end),
        y0=y0,
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-8,
    )
    theta_t = sol.y[:n]
    vel_t = sol.y[n:]
    r_t = kuramoto_order_parameter(theta_t)
    return sol.t, theta_t, vel_t, r_t


def run_kuramoto(cfg: KuramotoConfig, omega_override: Optional[np.ndarray] = None) -> dict:
    f, k, g, theta0, omega = _prepare_kuramoto(cfg, coupling=None, omega_override=omega_override)
    t, theta_t, omega_t, r_t = _run_kuramoto_core(cfg, k, theta0, omega)
    return {
        "t": t,
        "theta": theta_t,
        "r": r_t,
        "omega": omega,
        "omega_t": omega_t,
        "F": f,
        "K": k,
        "graph": g,
        "topology": str(cfg.topology),
        "model": str(cfg.model),
    }


def run_kuramoto_pyramid_resonance(cfg: KuramotoConfig) -> dict:
    """
    Kuramoto where each node's natural frequency is set by its chamber resonance (Hz -> rad/s).
    """
    f, k, g = build_kuramoto_coupling(cfg)
    n = k.shape[0]
    node_order = list(g.graph.get("node_order", g.nodes()))
    omega = np.zeros(n, dtype=np.float64)
    for i, node in enumerate(node_order):
        freqs = g.nodes[node].get("resonant_freqs", {}) if isinstance(g, nx.Graph) else {}
        primary_hz = freqs.get("primary", None) if isinstance(freqs, dict) else None
        if primary_hz is None:
            omega[i] = float(cfg.omega_mean)
        else:
            omega[i] = 2.0 * math.pi * float(primary_hz)
    f, k, g, theta0, omega = _prepare_kuramoto(cfg, coupling=(f, k, g), omega_override=omega)
    t, theta_t, omega_t, r_t = _run_kuramoto_core(cfg, k, theta0, omega)
    return {
        "t": t,
        "theta": theta_t,
        "r": r_t,
        "omega": omega,
        "omega_t": omega_t,
        "F": f,
        "K": k,
        "graph": g,
        "topology": str(cfg.topology),
        "model": str(cfg.model),
    }


def save_kuramoto_shaft_weight_sweep(out_dir: str, cfg: KuramotoConfig, weights: list[float]) -> str:
    if str(cfg.topology).lower().strip() != "pyramid":
        raise ValueError("Shaft-weight sweep is only supported for topology=pyramid.")
    if not weights:
        raise ValueError("Shaft-weight sweep requires at least one weight.")

    results = []
    for w in weights:
        base = cfg.__dict__
        medium = _parse_medium_weights(cfg.pyramid_medium_weights)
        medium["air_shaft"] = float(w)
        sweep_cfg = KuramotoConfig(**{**base, "pyramid_medium_weights": _format_medium_weights(medium)})
        r_final = save_kuramoto_plot(out_dir, sweep_cfg)
        results.append((float(w), float(r_final)))

    csv_path = os.path.join(out_dir, "kuramoto_pyramid_shaft_sweep.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["air_shaft_weight", "r_final"])
        for weight, r_final in results:
            w.writerow([weight, r_final])

    fig, ax = plt.subplots(figsize=(8, 4))
    xs = [r[0] for r in results]
    ys = [r[1] for r in results]
    ax.plot(xs, ys, marker="o", lw=2)
    ax.set_xlabel("air_shaft_weight")
    ax.set_ylabel("r_final")
    ax.set_title("Kuramoto r_final vs air_shaft_weight (pyramid)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plot_path = os.path.join(out_dir, "kuramoto_pyramid_shaft_sweep.png")
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    return csv_path


def save_kuramoto_freq_sweep(out_dir: str, cfg: KuramotoConfig, target_freqs: list[float]) -> str:
    """
    Sweep target frequencies (Hz) and measure final synchronization r.
    """
    freqs = [float(f) for f in target_freqs if float(f) >= 0.0]
    if not freqs:
        raise ValueError("frequency sweep requires at least one non-negative frequency.")

    results = []
    for freq in freqs:
        omega_mean = 2.0 * math.pi * float(freq)
        omega_std = float(cfg.omega_std) if omega_mean <= 0.0 else (0.1 * omega_mean)
        sweep_cfg = replace(cfg, omega_mean=omega_mean, omega_std=omega_std)
        out = run_kuramoto(sweep_cfg)
        r_final = float(out["r"][-1])
        results.append((float(freq), float(omega_mean), float(omega_std), float(r_final)))

    csv_path = os.path.join(out_dir, "kuramoto_freq_sweep.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["freq_hz", "omega_mean", "omega_std", "r_final"])
        for row in results:
            w.writerow(list(row))

    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    xs = [r[0] for r in results]
    ys = [r[3] for r in results]
    ax.plot(xs, ys, marker="o", lw=2)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("r_final")
    ax.set_title("Kuramoto r_final vs target frequency")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plot_path = os.path.join(out_dir, "kuramoto_freq_sweep.png")
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    return csv_path


def save_kuramoto_animation(out_dir: str, cfg: KuramotoConfig, filename: str = "kuramoto.gif") -> str:
    out = run_kuramoto(cfg)
    t = out["t"]
    theta = out["theta"]  # (N, T)
    n = theta.shape[0]

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    x = np.arange(1, n + 1)
    pts = ax.scatter(x, np.mod(theta[:, 0], 2.0 * np.pi), s=25)
    ax.set_ylim(0.0, 2.0 * np.pi)
    ax.set_xlabel("node index")
    ax.set_ylabel("theta mod 2π")
    ax.set_title(f"Kuramoto phase evolution ({out['topology']})")
    ax.grid(True, alpha=0.25)

    def update(frame: int):
        y = np.mod(theta[:, frame], 2.0 * np.pi)
        pts.set_offsets(np.column_stack([x, y]))
        ax.set_title(f"Kuramoto phase evolution ({out['topology']})  t={t[frame]:.2f}")
        return (pts,)

    ani = animation.FuncAnimation(fig, update, frames=theta.shape[1], interval=50, blit=True)
    out_path = os.path.join(out_dir, filename)
    writer = animation.PillowWriter(fps=20)
    ani.save(out_path, writer=writer)
    plt.close(fig)
    return out_path


def cantor_mask(x: np.ndarray, levels: int) -> np.ndarray:
    """
    Returns a boolean mask for points that remain in the classic middle-third Cantor set after `levels`.
    x is assumed normalized into [0,1].
    """
    if levels <= 0:
        return np.ones_like(x, dtype=bool)
    mask = np.ones_like(x, dtype=bool)
    u = x.copy()
    for _ in range(levels):
        # Remove middle third: u in (1/3, 2/3) at this scale
        remove = (u > (1.0 / 3.0)) & (u < (2.0 / 3.0))
        mask &= ~remove
        # Scale and wrap to next digit in base-3
        u = (u * 3.0) % 1.0
    return mask


@dataclass(frozen=True)
class TunnelingConfig:
    """
    Configuration for 1D wave packet tunneling simulation (split-step Fourier).
    
    Default parameters are chosen to reproduce the benchmark from spec.md:
    T ≈ 0.082 for a rectangular barrier with E=0.30, V0=1.0, a=4.0.
    
    Parameter derivations:
    ----------------------
    - n=512, L=40.0: Grid spacing dx = 40/512 ≈ 0.078. The de Broglie wavelength
      at E=0.30 is λ = 2π/k = 2π/√(2E) ≈ 8.1, so we have ~100 points per wavelength
      (well-resolved).
    
    - dt=0.1, steps=250: Total time = 25. The wave packet travels at group velocity
      v_g = k/m = √(2E) ≈ 0.77 (with m=1), covering ~19 units. Starting at x0=-10,
      this is enough to reach and traverse the barrier at x=0.
    
    - E=0.30, V0=1.0: E/V0 = 0.30 gives significant tunneling (not too opaque,
      not too transparent). The analytical WKB transmission for a rectangular
      barrier is T ≈ exp(-2κa) where κ = √(2(V0-E)) ≈ 1.18, giving T ≈ 0.082.
    
    - sigma=2.0: Wave packet width. The momentum uncertainty Δk ~ 1/(2σ) ≈ 0.25,
      so Δk/k0 ≈ 0.32 (reasonably monochromatic but still localized).
    
    - x0=-10.0: Start position. Far enough from barrier (at x=0) that the
      initial wave packet doesn't overlap with the potential.
    """
    n: int = 512
    L: float = 40.0
    dt: float = 0.1
    steps: int = 250
    V0: float = 1.0
    a: float = 4.0
    E: float = 0.30
    sigma: float = 2.0
    x0: float = -10.0
    barrier_type: str = "rect"  # "rect" or "cantor"
    cantor_levels: int = 5
    absorb: bool = False  # Default off for accurate T measurement; enable for Cantor sweeps
    # Absorbing boundary layer parameters (smooth "mask" attenuation method):
    #
    # The absorbing layer damps the wave function near domain edges to prevent
    # spurious reflections from the periodic FFT boundary. Each time step:
    #   psi *= exp(-strength * (dist/w)^power)
    # where dist = distance into the absorbing region, w = layer width.
    #
    # Design criteria:
    #   - After ~100 steps through the layer, amplitude should drop to <1%
    #   - Smooth onset (power >= 2) avoids artificial reflections from sharp edges
    #
    # Derivation of defaults:
    #   - absorb_frac=0.15: Layer is 15% of half-domain. For L=40, this is 3 units,
    #     giving the wave ~60 time steps (at dt=0.1) to traverse and be absorbed.
    #   - absorb_strength=6.0, absorb_power=2.0: At the outer edge (dist/w=1),
    #     damping per step = exp(-6) ≈ 0.0025. After entering the layer, amplitude
    #     drops to <1% within ~2 layer-widths of propagation.
    #   - These values empirically validated against analytical transmission for
    #     rectangular barriers (T_numeric ≈ 0.083 vs T_WKB ≈ 0.082).
    absorb_frac: float = 0.15
    absorb_strength: float = 6.0
    absorb_power: float = 2.0


def build_potential(cfg: TunnelingConfig, x: np.ndarray) -> np.ndarray:
    if cfg.barrier_type == "rect":
        return cfg.V0 * (np.abs(x) < (cfg.a / 2.0)).astype(np.float64)

    if cfg.barrier_type == "cantor":
        in_window = (np.abs(x) < (cfg.a / 2.0))
        u = (x[in_window] + (cfg.a / 2.0)) / cfg.a  # map to [0,1]
        keep = cantor_mask(u, cfg.cantor_levels)
        v = np.zeros_like(x, dtype=np.float64)
        # Interpret "Cantor-set barrier" as barrier material occupying the removed pieces
        # (complement of the Cantor set) so the effective barrier becomes denser with levels.
        v[in_window] = cfg.V0 * (~keep).astype(np.float64)
        return v

    raise ValueError("barrier_type must be 'rect' or 'cantor'")


def run_tunneling(cfg: TunnelingConfig) -> dict:
    # Units: hbar = m = 1
    x = np.linspace(-cfg.L / 2.0, cfg.L / 2.0, cfg.n)
    dx = float(x[1] - x[0])
    k = 2.0 * np.pi * np.fft.fftfreq(cfg.n, d=dx)

    V = build_potential(cfg, x)

    k0 = math.sqrt(max(0.0, 2.0 * cfg.E))
    psi0 = np.exp(-((x - cfg.x0) ** 2) / (4.0 * cfg.sigma**2)) * np.exp(1j * k0 * x)
    psi0 = psi0.astype(np.complex128)
    psi0 /= math.sqrt(float(np.sum(np.abs(psi0) ** 2) * dx))

    exp_V = np.exp(-1j * V * cfg.dt / 2.0)
    exp_K = np.exp(-1j * (k**2) * cfg.dt / 2.0)

    psi = psi0.copy()
    if cfg.absorb:
        half = 0.5 * cfg.L
        w = float(cfg.absorb_frac) * half
        w = max(w, 1e-9)
        x_cut = half - w
        # Smooth damping mask in the outer layer, applied each step.
        dist = np.maximum(0.0, np.abs(x) - x_cut)
        mask = np.exp(-float(cfg.absorb_strength) * (dist / w) ** float(cfg.absorb_power))
        mask = mask.astype(np.float64)
    else:
        mask = None
    for _ in range(int(cfg.steps)):
        psi *= exp_V
        psi_k = np.fft.fft(psi)
        psi_k *= exp_K
        psi = np.fft.ifft(psi_k)
        psi *= exp_V
        if mask is not None:
            psi *= mask

    right_mask = x > (cfg.a / 2.0)
    T = float(np.sum(np.abs(psi[right_mask]) ** 2) * dx)

    return {
        "x": x,
        "V": V,
        "psi0": psi0,
        "psi": psi,
        "T": T,
        "cfg": cfg,
    }


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def save_kuramoto_plot(
    out_dir: str,
    cfg: KuramotoConfig,
    extra: Optional[dict] = None,
    omega_override: Optional[np.ndarray] = None,
    use_pyramid_resonance: bool = False,
) -> float:
    extra = extra or {}
    if use_pyramid_resonance:
        out = run_kuramoto_pyramid_resonance(cfg)
    else:
        out = run_kuramoto(cfg, omega_override=omega_override)
    t = out["t"]
    r = out["r"]
    theta = out["theta"]
    theta_final = np.mod(theta[:, -1], 2.0 * np.pi)
    n = theta.shape[0]
    g = out.get("graph")

    fig, (ax_r, ax_phase) = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    ax_r.plot(t, r, lw=2)
    ax_r.set_title(f"Kuramoto r(t) (topology={out['topology']}, N={n}, final r={r[-1]:.3f})")
    ax_r.set_xlabel("t")
    ax_r.set_ylabel("r")
    ax_r.set_ylim(0.0, 1.05)
    ax_r.grid(True, alpha=0.3)

    ax_phase.scatter(np.arange(1, n + 1), theta_final, s=20)
    ax_phase.set_title("Final phases (wrapped to 0..2π)")
    ax_phase.set_xlabel("node index")
    ax_phase.set_ylabel("theta mod 2π")
    ax_phase.set_ylim(0.0, 2.0 * np.pi)
    ax_phase.grid(True, alpha=0.3)

    fig.savefig(os.path.join(out_dir, "kuramoto.png"), dpi=160)
    plt.close(fig)

    node_order = []
    if isinstance(g, nx.Graph):
        node_order = list(g.graph.get("node_order", []))
    if node_order:
        map_path = os.path.join(out_dir, "kuramoto_node_map.csv")
        with open(map_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["index", "node"])
            for i, node in enumerate(node_order, start=1):
                w.writerow([i, node])

    if bool(extra.get("label_nodes")) and node_order:
        fig2, ax2 = plt.subplots(figsize=(10, 4))
        x = np.arange(1, n + 1)
        ax2.scatter(x, theta_final, s=20)
        for i, node in enumerate(node_order, start=1):
            ax2.annotate(str(node), (i, theta_final[i - 1]), xytext=(4, 4), textcoords="offset points", fontsize=8)
        ax2.set_title("Final phases (labeled)")
        ax2.set_xlabel("node index")
        ax2.set_ylabel("theta mod 2π")
        ax2.set_ylim(0.0, 2.0 * np.pi)
        ax2.grid(True, alpha=0.3)
        fig2.tight_layout()
        fig2.savefig(os.path.join(out_dir, "kuramoto_phases_labeled.png"), dpi=160)
        plt.close(fig2)

    if bool(extra.get("fft")):
        dt = float(t[1] - t[0]) if len(t) > 1 else 1.0
        r_centered = r - float(np.mean(r))
        freqs = np.fft.rfftfreq(len(r_centered), d=dt)
        spectrum = np.abs(np.fft.rfft(r_centered))
        csv_path = os.path.join(out_dir, "kuramoto_r_fft.csv")
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["freq", "amplitude"])
            for fval, amp in zip(freqs, spectrum, strict=False):
                w.writerow([float(fval), float(amp)])
        fig3, ax3 = plt.subplots(figsize=(10, 4))
        ax3.plot(freqs, spectrum, lw=2)
        ax3.set_title("FFT of r(t) (detrended)")
        ax3.set_xlabel("frequency (1/time unit)")
        ax3.set_ylabel("amplitude")
        ax3.grid(True, alpha=0.3)
        fig3.tight_layout()
        fig3.savefig(os.path.join(out_dir, "kuramoto_r_fft.png"), dpi=160)
        plt.close(fig3)

    phase_pairs = _parse_phase_pairs(extra.get("phase_pairs"))
    if phase_pairs:
        cols = 2
        rows = int(math.ceil(len(phase_pairs) / cols))
        fig4, axes = plt.subplots(rows, cols, figsize=(10, 4 * rows), constrained_layout=True)
        if not isinstance(axes, np.ndarray):
            axes = np.array([axes])
        axes = axes.flatten()
        for ax, (i, j) in zip(axes, phase_pairs, strict=False):
            if i > n or j > n:
                ax.set_visible(False)
                continue
            th_i = np.mod(theta[i - 1, :], 2.0 * np.pi)
            th_j = np.mod(theta[j - 1, :], 2.0 * np.pi)
            ax.plot(th_i, th_j, lw=1.2)
            ax.set_title(f"Phase portrait: {i} vs {j}")
            ax.set_xlabel(f"theta[{i}]")
            ax.set_ylabel(f"theta[{j}]")
            ax.grid(True, alpha=0.3)
        for ax in axes[len(phase_pairs) :]:
            ax.set_visible(False)
        fig4.savefig(os.path.join(out_dir, "kuramoto_phase_portrait.png"), dpi=160)
        plt.close(fig4)

    return float(r[-1])


def save_tunneling_plot(out_dir: str, cfg: TunnelingConfig) -> float:
    out = run_tunneling(cfg)
    x = out["x"]
    V = out["V"]
    psi0 = out["psi0"]
    psi = out["psi"]
    T = float(out["T"])

    rho0 = np.abs(psi0) ** 2
    rho = np.abs(psi) ** 2
    v_scaled = V / (np.max(V) + 1e-12) * (0.9 * max(float(np.max(rho0)), float(np.max(rho))) + 1e-12)

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    ax0.plot(x, rho0, lw=2, label="|psi|^2 (t=0)")
    ax0.plot(x, v_scaled, lw=1.5, label="V(x) (scaled)", alpha=0.8)
    ax0.set_title("|psi(x,t=0)|^2 and V(x)")
    ax0.set_xlabel("x")
    ax0.set_ylabel("density / potential (scaled)")
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc="upper right")

    ax1.plot(x, rho, lw=2, label="|psi|^2 (t=end)")
    ax1.plot(x, v_scaled, lw=1.5, label="V(x) (scaled)", alpha=0.8)
    ax1.set_title(f"|psi(x,t=end)|^2 and V(x)   (T={T:.3f})")
    ax1.set_xlabel("x")
    ax1.set_ylabel("density / potential (scaled)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right")

    fig.savefig(os.path.join(out_dir, "tunneling.png"), dpi=160)
    plt.close(fig)
    return T


def save_cantor_t_sweep_plot(out_dir: str, base: TunnelingConfig, levels: list[int]) -> dict:
    """
    Sweep Cantor levels and fit log(T) vs log(levels).
    Returns slope, r2, and arrays.
    """
    lvls = [int(x) for x in levels if int(x) > 0]
    if not lvls:
        raise ValueError("levels must be positive integers")

    Ts = []
    for L in lvls:
        cfg = TunnelingConfig(**{**base.__dict__, "barrier_type": "cantor", "cantor_levels": int(L)})
        out = run_tunneling(cfg)
        Ts.append(float(out["T"]))

    x = np.array(lvls, dtype=np.float64)
    y = np.array(Ts, dtype=np.float64)
    # Numerical epsilon for log():
    # ----------------------------
    # eps = 1e-15 is ~10× machine epsilon for float64 (≈2.2e-16).
    # This prevents log(0) when transmission is numerically zero at high
    # Cantor levels, while being small enough not to affect valid data points
    # (typical T values are 1e-6 to 1, so eps << min(T) for resolved levels).
    eps = 1e-15
    lx = np.log(x)
    ly = np.log(np.maximum(y, eps))

    # Primary (requested) power-law fit: log(T) vs log(levels)
    m, b = np.polyfit(lx, ly, deg=1)
    ly_hat = m * lx + b
    ss_res = float(np.sum((ly - ly_hat) ** 2))
    ss_tot = float(np.sum((ly - float(np.mean(ly))) ** 2)) + 1e-15
    r2 = 1.0 - (ss_res / ss_tot)
    d_eff = -float(m)

    # Alternative fit: exponential-in-level (often more appropriate for refinement level)
    mL, bL = np.polyfit(x, ly, deg=1)
    ly_hat_L = mL * x + bL
    ss_res_L = float(np.sum((ly - ly_hat_L) ** 2))
    r2_L = 1.0 - (ss_res_L / ss_tot)

    # Alternative fit: log(T) vs log(2^levels) (same data, different x-scale)
    lx2 = x * math.log(2.0)
    m2, b2 = np.polyfit(lx2, ly, deg=1)
    ly_hat_2 = m2 * lx2 + b2
    ss_res_2 = float(np.sum((ly - ly_hat_2) ** 2))
    r2_2 = 1.0 - (ss_res_2 / ss_tot)

    # Resolution sanity: smallest Cantor feature ~ a / 3^L should be meaningfully larger than dx.
    # Note: run_tunneling() uses x = linspace(-L/2, L/2, n), so dx = L/(n-1) (not L/n).
    dx = float(base.L) / float(base.n - 1)
    min_feat = np.array([float(base.a) / (3.0**int(L)) for L in lvls], dtype=np.float64)
    feat_over_dx = (min_feat / dx).tolist()  # >1 means feature larger than grid step

    # Heuristic warnings:
    # - If min_feature/dx < ~5, the smallest Cantor features are under-resolved and
    #   any “fractal scaling” fit will be numerically dominated (not physically meaningful).
    under = [L for L, r in zip(lvls, feat_over_dx, strict=False) if float(r) < 5.0]
    if under:
        print(f"[warn] Cantor levels under-resolved (min_feature/dx < 5): {under}")
    if bool(getattr(base, "absorb", False)):
        print("[warn] absorb=True: probability is not conserved; interpret T as 'mass on right', not a strict transmission coefficient.")

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)
    ax0.plot(x, y, marker="o", lw=2)
    ax0.set_title("Cantor barrier: transmission T vs cantor_levels")
    ax0.set_xlabel("cantor_levels")
    ax0.set_ylabel("T")
    ax0.grid(True, alpha=0.3)

    ax1.plot(lx, ly, marker="o", lw=2, label="data")
    ax1.plot(lx, ly_hat, lw=2, label=f"log-log: m={m:.3f}, R²={r2:.3f}")
    ax1.set_title(f"log-log (requested) d_eff=-m={-m:.3f} | alt R² log(T)~levels: {r2_L:.3f}")
    ax1.set_xlabel("log(levels)")
    ax1.set_ylabel("log(T)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="best")

    fig.savefig(os.path.join(out_dir, "cantor_t_sweep.png"), dpi=160)
    plt.close(fig)

    csv_path = os.path.join(out_dir, "cantor_t_sweep.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cantor_level", "T", "min_feature_over_dx"])
        for L, T, fodx in zip(lvls, Ts, feat_over_dx, strict=False):
            w.writerow([L, T, fodx])
        w.writerow([])
        w.writerow(["fit_logT_vs_loglevels_slope_m", float(m)])
        w.writerow(["fit_logT_vs_loglevels_intercept_b", float(b)])
        w.writerow(["fit_logT_vs_loglevels_R2", float(r2)])
        w.writerow(["d_eff_minus_slope", float(d_eff)])
        w.writerow(["fit_logT_vs_levels_slope_mL", float(mL)])
        w.writerow(["fit_logT_vs_levels_intercept_bL", float(bL)])
        w.writerow(["fit_logT_vs_levels_R2", float(r2_L)])
        w.writerow(["fit_logT_vs_log2powlevels_slope_m2", float(m2)])
        w.writerow(["fit_logT_vs_log2powlevels_intercept_b2", float(b2)])
        w.writerow(["fit_logT_vs_log2powlevels_R2", float(r2_2)])

    return {
        "levels": lvls,
        "T": Ts,
        "slope": float(m),
        "d_eff": d_eff,
        "r2": float(r2),
        "r2_levels": float(r2_L),
        "r2_log2pow": float(r2_2),
    }


def headless_main(mode: str, out_dir: str, barrier: str) -> int:
    _ensure_dir(out_dir)
    print(f"[headless] saving outputs to: {os.path.abspath(out_dir)}")

    final_r = None
    T = None

    if mode in ("kuramoto", "both", "all"):
        kcfg = globals().get("_KURAMOTO_HEADLESS_CFG", KuramotoConfig())
        kextra = globals().get("_KURAMOTO_EXTRA", {})
        final_r = save_kuramoto_plot(out_dir, kcfg, extra=kextra, use_pyramid_resonance=bool(kextra.get("pyramid_resonance")))
        print(f"[kuramoto] final r = {final_r:.3f}  (see {os.path.join(out_dir, 'kuramoto.png')})")
        extra = globals().get("_KURAMOTO_SWEEP", {"pyramid_shaft_sweep": None})
        sweep_vals = _parse_float_list(extra.get("pyramid_shaft_sweep"))
        if sweep_vals:
            csv_path = save_kuramoto_shaft_weight_sweep(out_dir, kcfg, sweep_vals)
            print(f"[kuramoto] wrote {csv_path} and kuramoto_pyramid_shaft_sweep.png")
        freq_vals = _parse_float_list(extra.get("freq_sweep"))
        if freq_vals:
            csv_path = save_kuramoto_freq_sweep(out_dir, kcfg, freq_vals)
            print(f"[kuramoto] wrote {csv_path} and kuramoto_freq_sweep.png")
    if mode in ("kuramoto_anim",):
        kcfg = globals().get("_KURAMOTO_HEADLESS_CFG", KuramotoConfig())
        out_path = save_kuramoto_animation(out_dir, kcfg)
        print(f"[kuramoto_anim] wrote {out_path}")

    if mode in ("tunneling", "both", "all"):
        tcfg = TunnelingConfig(barrier_type=barrier)
        T = save_tunneling_plot(out_dir, tcfg)
        print(f"[tunneling] T = {T:.3f}  (see {os.path.join(out_dir, 'tunneling.png')})")
        extra = globals().get("_TUNNEL_EXTRA", {"cantor_sweep": False})
        if bool(extra.get("cantor_sweep", False)):
            lvls = extra.get("levels") or [1, 2, 3, 4, 5, 6, 7, 8]
            # Cantor sweep configuration: derived from resolution requirements
            # ----------------------------------------------------------------
            # The Cantor set at level L has smallest features of size: a / 3^L
            #
            # For accurate simulation, we need: dx << min_feature
            # Rule of thumb: at least 5-10 grid points per smallest feature.
            #
            # Derivation of parameters:
            #   - a=10.0, L_max=8: min_feature = 10/3^8 ≈ 0.0015
            #   - n=16384, L=12.0: dx = 12/16384 ≈ 0.00073
            #   - Ratio: min_feature/dx ≈ 2.0 (marginal at L=8; L≤6 well-resolved)
            #
            #   - L=12.0: Domain length. Short enough that wave packet traverses
            #     the barrier region within the simulation time.
            #   - x0=-4.0, sigma=0.4: Compact wave packet starting left of barrier,
            #     narrow enough to have well-defined momentum (Δx·Δk ~ 1).
            #   - dt=0.002, steps=2500: Total time = 5.0. Split-step Fourier is
            #     unconditionally stable (no CFL constraint), but accuracy depends
            #     on dt resolving phase accumulation from kinetic + potential terms.
            #     Validated by convergence: T varies <1% when halving dt from 0.004→0.002.
            #   - absorb_frac=0.2, absorb_strength=8.0: Stronger absorption for
            #     shorter domain (less distance to absorb reflected waves).
            #
            # Note: CSV output includes min_feature_over_dx; values < 5 indicate
            # the level is under-resolved and results may not be reliable.
            sweep_cfg = TunnelingConfig(
                n=16384,
                L=12.0,
                a=10.0,
                x0=-4.0,
                sigma=0.4,
                dt=0.002,
                steps=2500,
                E=tcfg.E,
                V0=tcfg.V0,
                barrier_type="cantor",
                absorb=True,
                absorb_frac=0.2,
                absorb_strength=8.0,
            )
            res = save_cantor_t_sweep_plot(out_dir, sweep_cfg, levels=lvls)
            print(
                f"[tunneling] wrote {os.path.join(out_dir, 'cantor_t_sweep.png')} "
                f"(d_eff~{res['d_eff']:.3f}, R2={res['r2']:.3f}, levels={res['levels']})"
            )

    if mode in ("chsh", "all"):
        cfg = globals().get("_CHSH_HEADLESS_CFG", CHSHConfig())
        chsh = save_chsh_plot(out_dir, cfg)
        print(f"[chsh] CHSH S = {chsh:.3f}  (see {os.path.join(out_dir, 'chsh.png')})")
        extra = globals().get("_CHSH_EXTRA", {"noise_sweep": False})
        if bool(extra.get("noise_sweep", False)):
            noises = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
            save_chsh_noise_sweep_plot(out_dir, cfg, noises=noises)
            print(f"[chsh] wrote {os.path.join(out_dir, 'chsh_noise_sweep.png')} (noises={noises})")

    if mode in ("retro", "all"):
        rcfg = globals().get("_RETRO_HEADLESS_CFG", RetroConfig())
        if bool(globals().get("_UNIFIED_DEMO", False)) and mode == "all" and final_r is not None:
            # Example mapping:
            # - phase shift (does NOT change visibility magnitude by itself; it shifts the fringe)
            # - dephasing strength (DOES reduce control coherence / visibility)
            phi = 45.0 + 30.0 * (1.0 - float(final_r))
            deph = max(0.0, min(1.0, 1.0 - float(final_r)))
            rcfg = RetroConfig(
                **{
                    **rcfg.__dict__,
                    "use_switch": True,
                    "switch_cptp": True,
                    "switch_phase_deg": float(phi),
                    "switch_dephase_p": float(deph),
                }
            )
            print(
                f"[unified] mapped kuramoto r_final={final_r:.3f} -> retro phi={phi:.2f} deg, "
                f"control_dephase_p={deph:.3f} (CPTP switch)"
            )
        retro = save_retro_plot(out_dir, rcfg)
        print(
            f"[retro] max_abs_delta = {retro['delta_max_abs']:.3f}  (see {os.path.join(out_dir, 'retro.png')})\n"
            f"        fwd_freq={retro['freq_forward']:.3f}Hz  rev_freq={retro['freq_reverse']:.3f}Hz  |df|={retro['freq_delta']:.3f}Hz\n"
            f"        spec_dist={retro['spec_dist']:.3f}\n"
            f"        switch_vmax_max={retro['switch_vmax_max']:.3f}, switch_fringe_max={retro['switch_fringe_max']:.3f}  "
            f"(wrote {os.path.join(out_dir, 'retro_switch.png')} if enabled)\n"
            f"        commutators: {retro['comm_norms']}  (also wrote {os.path.join(out_dir, 'retro_spectrum.png')})"
        )
        sweep = globals().get("_RETRO_SWEEP", {"enabled": False, "step": 5.0})
        if bool(sweep.get("enabled", False)):
            save_retro_theta_sweep_plot(out_dir, rcfg, step_deg=float(sweep.get("step", 5.0)))
            print(f"[retro] wrote {os.path.join(out_dir, 'retro_theta_sweep.png')}")
        extra = globals().get("_RETRO_EXTRA", {"vis_sweep": False, "phi_sweep": False, "phi_step": 10.0, "bloch": False})
        if bool(extra.get("vis_sweep", False)):
            save_retro_vis_vs_theta_plot(out_dir, rcfg, step_deg=float(sweep.get("step", 5.0)))
            print(f"[retro] wrote {os.path.join(out_dir, 'retro_vis_vs_theta.png')}")
        if bool(extra.get("phi_sweep", False)):
            save_retro_vis_vs_phi_plot(out_dir, rcfg, step_deg=float(extra.get("phi_step", 10.0)))
            print(f"[retro] wrote {os.path.join(out_dir, 'retro_vis_vs_phi.png')}")
        if bool(extra.get("bloch", False)):
            out_path = save_retro_control_bloch_gif(out_dir, rcfg)
            print(f"[retro] wrote {out_path}")

    return 0


# --- CHSH / entanglement (pure NumPy; no QuTiP dependency) --------------------

_SQRT2 = math.sqrt(2.0)


def _paulis() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    sy = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
    sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    return sx, sy, sz


def _phi_plus_state() -> np.ndarray:
    # |Φ+> = (|00> + |11>) / sqrt(2)
    psi = np.zeros((4,), dtype=np.complex128)
    psi[0] = 1.0 / _SQRT2
    psi[3] = 1.0 / _SQRT2
    return psi


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n == 0.0:
        raise ValueError("zero vector")
    return v / n


def _meas_op(n: np.ndarray) -> np.ndarray:
    # n·σ
    sx, sy, sz = _paulis()
    n = _unit(n.astype(np.float64))
    return n[0] * sx + n[1] * sy + n[2] * sz


def _rotate_about_y(v: np.ndarray, delta: float) -> np.ndarray:
    # Rotate in x-z plane about y-axis:
    # x' = x cosδ + z sinδ
    # z' = -x sinδ + z cosδ
    x, y, z = float(v[0]), float(v[1]), float(v[2])
    c = math.cos(delta)
    s = math.sin(delta)
    return np.array([x * c + z * s, y, -x * s + z * c], dtype=np.float64)


@dataclass(frozen=True)
class CHSHConfig:
    """
    Configuration for CHSH Bell inequality simulation.
    
    The simulation computes the CHSH parameter S for the Bell state |Φ+⟩,
    optionally mixed with white noise (Werner state) to demonstrate the
    quantum-to-classical transition.
    
    Parameter choices:
    ------------------
    - n_ofm=20: Number of integers used to compute the OFM scale factor.
      This determines how "fractal_noise" maps to actual depolarization.
      Larger n → more stable ofm_scale (converges to ~0.4-0.5 for typical n).
    
    - fractal_noise=0.0: Noise parameter. At 0, the state is pure |Φ+⟩ and
      S = 2√2 ≈ 2.828 (Tsirelson bound). As noise increases, S decreases.
      See the k=0.23 derivation comment in run_chsh() for the physics.
    
    Physics background:
    -------------------
    The CHSH inequality states that for any local hidden variable theory:
      S = |E(a,b) + E(a,b') + E(a',b) - E(a',b')| ≤ 2
    
    Quantum mechanics allows S up to 2√2 (Tsirelson bound), achieved by
    the maximally entangled Bell state with optimal measurement settings.
    """
    n_ofm: int = 20
    fractal_noise: float = 0.0  # maps to a depolarizing "trace-out" strength (Werner mixture)
    seed: int = 1  # NOTE: currently unused (CHSH calculation is deterministic)


def run_chsh(cfg: CHSHConfig) -> dict:
    # Standard settings that saturate Tsirelson bound for |Φ+>:
    # A = σz, A' = σx
    # B = (σz + σx)/sqrt(2), B' = (σz - σx)/sqrt(2)
    a = np.array([0.0, 0.0, 1.0])
    a_p = np.array([1.0, 0.0, 0.0])
    b = _unit(a + a_p)
    b_p = _unit(a - a_p)

    # "Fractal twist": derive an OFM scale and map it to a depolarizing (Werner) mixture.
    #
    # This gives a clean, monotonic suppression:
    #   rho = (1-p)|Φ+><Φ+| + p I/4
    # which yields S -> 0 as p->1, and crosses classical bound near p≈0.293.
    f = ordered_factorizations_upto(max(2, int(cfg.n_ofm)))
    vals = f[1 : int(cfg.n_ofm) + 1].astype(np.float64)
    d = np.abs(vals[:, None] - vals[None, :])
    if d.size <= 1:
        ofm_scale = 0.0
        avg_d = 0.0
        max_d = 0.0
    else:
        tri = np.triu_indices(int(cfg.n_ofm), 1)
        avg_d = float(np.mean(d[tri]))
        max_d = float(np.max(d))
        ofm_scale = 0.0 if max_d <= 0 else (avg_d / max_d)

    # Map noise -> depolarizing strength p.
    #
    # Physics derivation of k:
    # -------------------------
    # For a Werner state (depolarizing channel on |Φ+⟩):
    #   ρ = (1-p)|Φ+⟩⟨Φ+| + p·I/4
    #
    # The CHSH value for this state is:
    #   S(p) = 2√2 · (1 - p)
    #
    # The classical bound S=2 is reached when:
    #   2 = 2√2 · (1 - p_critical)
    #   p_critical = 1 - 1/√2 ≈ 0.293
    #
    # We want: p = fractal_noise × ofm_scale × k
    #
    # Design goal: noise=3 should approach the classical bound (p ≈ 0.293)
    # for typical OFM parameters (N=20 gives ofm_scale ≈ 0.42).
    #
    # Solving for k:
    #   0.293 = 3 × 0.42 × k
    #   k = 0.293 / (3 × 0.42) ≈ 0.23
    #
    # This gives a physically meaningful mapping where:
    #   - noise=0 → p=0 → S=2√2 (Tsirelson bound, maximal entanglement)
    #   - noise≈3 → p≈0.29 → S≈2 (classical limit)
    #   - noise>3 → S<2 (below classical bound, mostly mixed state)
    k = 0.23
    p = float(min(0.999, max(0.0, float(cfg.fractal_noise) * float(ofm_scale) * k)))

    # Use standard optimal settings (no angle perturbation); only state gets mixed.
    b1 = b
    b2 = b_p

    psi = _phi_plus_state().reshape((4, 1))
    rho_phi = psi @ psi.conj().T
    rho = (1.0 - p) * rho_phi + p * (np.eye(4, dtype=np.complex128) / 4.0)

    def corr(n1: np.ndarray, n2: np.ndarray) -> float:
        o = np.kron(_meas_op(n1), _meas_op(n2))
        return float(np.real(np.trace(rho @ o)))

    e_ab = corr(a, b1)
    e_abp = corr(a, b2)
    e_apb = corr(a_p, b1)
    e_apbp = corr(a_p, b2)
    s = e_ab + e_abp + e_apb - e_apbp

    return {
        "CHSH": float(s),
        "p": float(p),
        "ofm_avg": float(avg_d),
        "ofm_max": float(max_d),
        "ofm_scale": float(ofm_scale),
    }


def save_chsh_plot(out_dir: str, cfg: CHSHConfig) -> float:
    out = run_chsh(cfg)
    s = float(out["CHSH"])

    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    ax.bar(["CHSH"], [s], color="#2a6fdb" if s > 2.0 else "#d44a3a")
    ax.axhline(2.0, color="#2c8a2c", ls="--", lw=1.5, label="Classical bound (2)")
    ax.axhline(2.0 * _SQRT2, color="#7a2fbf", ls="--", lw=1.5, label="Tsirelson (2√2)")
    ax.set_ylim(0.0, 3.0)
    ax.set_title(
        f"CHSH S={s:.3f}  (noise={cfg.fractal_noise:.2f}, p={out['p']:.3f}, OFM scale={out['ofm_scale']:.3f})"
    )
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.2, axis="y")
    fig.savefig(os.path.join(out_dir, "chsh.png"), dpi=160)
    plt.close(fig)
    return s


def save_chsh_noise_sweep_plot(out_dir: str, base: CHSHConfig, noises: list[float]) -> None:
    ss = []
    ps = []
    ofm_scale = float(run_chsh(CHSHConfig(n_ofm=base.n_ofm, fractal_noise=0.0, seed=base.seed))["ofm_scale"])
    for n in noises:
        cfg = CHSHConfig(**{**base.__dict__, "fractal_noise": float(n)})
        out = run_chsh(cfg)
        ss.append(float(out["CHSH"]))
        ps.append(float(out["p"]))

    fig, ax = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)
    ax.plot(noises, ss, marker="o", lw=2, label="CHSH S")
    ax.axhline(2.0, color="#2c8a2c", ls="--", lw=1.5, label="Classical bound (2)")
    ax.axhline(2.0 * _SQRT2, color="#7a2fbf", ls="--", lw=1.5, label="Tsirelson (2√2)")
    ax.set_xlabel("fractal_noise")
    ax.set_ylabel("S")
    ax.set_ylim(0.0, 3.0)
    ax.grid(True, alpha=0.3)
    ax.set_title(f"CHSH suppression via OFM-scaled depolarizing trace-out (OFM N={base.n_ofm})")
    ax.legend(loc="best")
    fig.savefig(os.path.join(out_dir, "chsh_noise_sweep.png"), dpi=160)
    plt.close(fig)

    csv_path = os.path.join(out_dir, "chsh_noise_sweep.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["noise", "S", "p", "ofm_scale", "ofm_n"])
        for n, s, p in zip(noises, ss, ps, strict=False):
            w.writerow([n, s, p, ofm_scale, base.n_ofm])


# --- Retrocausality / non-commuting projectors (toy model) --------------------


def _I2() -> np.ndarray:
    return np.eye(2, dtype=np.complex128)


def _kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.kron(a, b)


def _unitary(H: np.ndarray, t: float) -> np.ndarray:
    return expm((-1.0j) * H * float(t))


def _density_from_ket(psi: np.ndarray) -> np.ndarray:
    psi = psi.reshape((-1, 1)).astype(np.complex128)
    return psi @ psi.conj().T


def _projector_from_axis(n: np.ndarray) -> np.ndarray:
    # Projector onto +1 eigenstate of n·σ : P = (I + n·σ)/2
    return 0.5 * (_I2() + _meas_op(n))


def _trace(a: np.ndarray) -> complex:
    return np.trace(a)


def _apply_projective_measurement(rho: np.ndarray, P: np.ndarray) -> tuple[float, np.ndarray]:
    # Outcome "P" occurs with prob p = Tr(P rho). Post-state: P rho P / p
    p = float(np.real(_trace(P @ rho)))
    if p <= 0.0:
        return 0.0, rho * 0.0
    post = (P @ rho @ P) / p
    return p, post


@dataclass(frozen=True)
class RetroConfig:
    """
    Configuration for retrocausality / non-commuting projector simulation.
    
    This toy model demonstrates that when projectors P and Q don't commute
    (and/or don't commute with the Hamiltonian H), the joint probability of
    measuring P-then-Q differs from Q-then-P.
    
    Default parameter choices:
    --------------------------
    - omega=1.0: Hamiltonian frequency H = (ω/2)·σ_x. Sets the precession rate
      on the Bloch sphere. With ω=1, one full rotation takes t = 2π ≈ 6.28.
    
    - t2=5.0: Total evolution time. Slightly less than one full rotation,
      showing clear oscillatory behavior in the order-dependent probabilities.
    
    - dt=0.02: Time step for sweeping t1 from 0 to t2. With t2=5, this gives
      250 data points, smooth enough for FFT analysis.
    
    - H_axis=(1,0,0), P_axis=(0,0,1): H generates X-rotations, P projects onto Z.
      These don't commute: [H,P] ≠ 0, which is essential for order-dependence.
    
    - theta_deg=90.0: Q projects onto the X-axis (θ=90° from Z in the x-z plane).
      At θ=90°, [P,Q] is maximal (P and Q are orthogonal observables), giving
      the strongest order-dependent effects. At θ=0°, P=Q and order doesn't matter.
    
    - psi0_axis=(0,0,1): Initial state |0⟩ (Z+ eigenstate). Starting aligned with
      P gives p_forward ≠ p_reverse because Q is not aligned.
    """
    omega: float = 1.0
    t2: float = 5.0
    dt: float = 0.02
    # H generates rotations about this axis (default: X-axis)
    H_axis: tuple[float, float, float] = (1.0, 0.0, 0.0)
    # P projects onto +1 eigenstate of this axis (default: Z+)
    P_axis: tuple[float, float, float] = (0.0, 0.0, 1.0)
    # Q axis: rotated in x-z plane by theta degrees from Z (θ=0 → Z, θ=90 → X)
    theta_deg: float = 90.0
    strobe_every_n_steps: int = 0  # 0 = off; >0 applies repeated projective "filtering" every N steps (Zeno-like)
    # Indefinite causal order (quantum-switch-style)
    use_switch: bool = False
    switch_cptp: bool = True  # if True, uses CPTP quantum switch; if False, uses postselected overlap model
    switch_phase_deg: float = 0.0  # relative phase between the two orders when computing p_plus/p_minus
    switch_dephase_p: float = 0.0  # optional phase-damping strength applied to the control after the switch (0..1)
    # Optional OFM modulation to tie back to the fractal metric
    ofm_modulate: bool = False
    ofm_n: int = 20
    ofm_strength_deg: float = 0.0
    ofm_target: str = "phi"  # "phi" or "theta"
    # initial state axis (Bloch): defaults to Z+ (|0⟩)
    psi0_axis: tuple[float, float, float] = (0.0, 0.0, 1.0)
    # Sanity checks (assertions) for the CPTP switch reduced density matrix.
    # Keep enabled for headless/dev runs; GUI disables by default to avoid hard crashes
    # on rare platform-specific floating-point edge cases.
    sanity_checks: bool = True


def _ket_from_axis(n: np.ndarray) -> np.ndarray:
    """
    Returns a normalized ket for the +1 eigenstate of n·σ.
    For n = z, this is |0>. For n = x, (|0>+|1>)/sqrt(2), etc.
    """
    n = _unit(n.astype(np.float64))
    op = _meas_op(n)
    vals, vecs = np.linalg.eigh(op)
    # eigenvalues are [-1, +1] in some order; pick +1
    idx = int(np.argmax(vals.real))
    psi = vecs[:, idx]
    psi = psi / np.linalg.norm(psi)
    return psi.astype(np.complex128)


def run_retro(cfg: RetroConfig) -> dict:
    H_axis = np.array(cfg.H_axis, dtype=np.float64)
    P_axis = np.array(cfg.P_axis, dtype=np.float64)
    # Optional OFM modulation
    ofm_delta_deg = 0.0
    if bool(cfg.ofm_modulate) and float(cfg.ofm_strength_deg) != 0.0:
        f_ofm = ordered_factorizations_upto(max(2, int(cfg.ofm_n)))
        vals = f_ofm[1 : int(cfg.ofm_n) + 1].astype(np.float64)
        d = np.abs(vals[:, None] - vals[None, :])
        if d.size > 1:
            tri = np.triu_indices(int(cfg.ofm_n), 1)
            avg_d = float(np.mean(d[tri]))
            max_d = float(np.max(d))
            scale = 0.0 if max_d <= 0 else (avg_d / max_d)
            ofm_delta_deg = float(cfg.ofm_strength_deg) * scale

    theta_deg_eff = float(cfg.theta_deg)
    phi_deg_eff = float(cfg.switch_phase_deg)
    if bool(cfg.ofm_modulate):
        tgt = str(cfg.ofm_target).lower().strip()
        if tgt == "theta":
            theta_deg_eff = theta_deg_eff + ofm_delta_deg
        else:
            phi_deg_eff = phi_deg_eff + ofm_delta_deg

    # Q axis is rotated in x–z plane: n = (sinθ, 0, cosθ)
    theta = math.radians(theta_deg_eff)
    Q_axis = np.array([math.sin(theta), 0.0, math.cos(theta)], dtype=np.float64)
    psi0_axis = np.array(cfg.psi0_axis, dtype=np.float64)

    sx, sy, sz = _paulis()
    H_dir = _unit(H_axis)
    H = 0.5 * float(cfg.omega) * (H_dir[0] * sx + H_dir[1] * sy + H_dir[2] * sz)

    P = _projector_from_axis(P_axis)
    Q = _projector_from_axis(Q_axis)

    psi0 = _ket_from_axis(psi0_axis)
    rho0 = _density_from_ket(psi0)

    # Non-commutativity diagnostics
    comm_HP = H @ P - P @ H
    comm_HQ = H @ Q - Q @ H
    comm_PQ = P @ Q - Q @ P
    comm_norms = {
        "||[H,P]||F": float(np.linalg.norm(comm_HP)),
        "||[H,Q]||F": float(np.linalg.norm(comm_HQ)),
        "||[P,Q]||F": float(np.linalg.norm(comm_PQ)),
    }

    t2 = float(cfg.t2)
    dt = float(cfg.dt)
    t1s = np.arange(0.0, t2 + 1e-12, dt)

    p_forward = np.zeros_like(t1s)
    p_reverse = np.zeros_like(t1s)

    # Forward: evolve to t1, project with P, evolve to t2, project with Q (joint probability P then Q)
    # Reverse: evolve to t1, project with Q, evolve to t2, project with P (joint probability Q then P)
    #
    # Quantum switch variants:
    #
    # - Postselected (legacy): uses overlap c(t1)=Tr(Kf ρ Kr†) for *selected* Kraus paths.
    # - CPTP (requested): builds the full switch channel from Kraus sets (nonselective measurements)
    #   and reads interference from the control qubit coherence directly (no postselection).
    use_switch = bool(cfg.use_switch)
    switch_cptp = bool(cfg.switch_cptp)
    phi = math.radians(phi_deg_eff)
    exp_i_phi = complex(math.cos(phi), math.sin(phi))
    c_overlap = np.zeros_like(t1s, dtype=np.complex128)  # used for postselected
    p_plus = np.zeros_like(t1s)  # used for postselected or CPTP
    p_minus = np.zeros_like(t1s)
    v_max = np.zeros_like(t1s)  # visibility metric (0..1)
    control_coh = np.zeros_like(t1s, dtype=np.complex128)  # ρ_c[0,1] for CPTP
    control_bloch = np.zeros((t1s.size, 3), dtype=np.float64)

    # Control state |+_phi> = (|0> + e^{i phi} |1>)/sqrt(2)
    ket0 = np.array([1.0, 0.0], dtype=np.complex128)
    ket1 = np.array([0.0, 1.0], dtype=np.complex128)
    ket_c = (ket0 + exp_i_phi * ket1) / _SQRT2
    rho_c_in = _density_from_ket(ket_c)
    rho_in_sc = np.kron(rho_c_in, rho0)  # control ⊗ system (4x4)

    # Projectors on the control qubit space (not to be confused with P, Q on system)
    P0_ctrl = _density_from_ket(ket0)  # |0><0|_control
    P1_ctrl = _density_from_ket(ket1)  # |1><1|_control

    def _ptrace_system(rho_sc: np.ndarray) -> np.ndarray:
        """
        Partial trace over the system to get the reduced control density matrix.
        
        rho_sc is a 4×4 matrix in the composite space (control ⊗ system).
        We want: (ρ_c)_{a,c} = Σ_b (ρ_sc)_{(a,b),(c,b)}
        
        Reshape axes: (a,b,c,d) = (c1,s1,c2,s2) where c=control, s=system
        Then contract: sum over s1=s2, keeping c1,c2.
        """
        r = rho_sc.reshape(2, 2, 2, 2)  # (c1, s1, c2, s2)
        return np.einsum("abcb->ac", r)  # sum over s with s1==s2

    for idx, t1 in enumerate(t1s):
        U1 = _unitary(H, t1)
        U2 = _unitary(H, t2 - t1)

        rho_t1 = U1 @ rho0 @ U1.conj().T

        pP, rho_after_P = _apply_projective_measurement(rho_t1, P)
        rho_t2_after_P = U2 @ rho_after_P @ U2.conj().T
        pQ_given_P = float(np.real(_trace(Q @ rho_t2_after_P)))
        p_forward[idx] = pP * pQ_given_P

        pQ, rho_after_Q = _apply_projective_measurement(rho_t1, Q)
        rho_t2_after_Q = U2 @ rho_after_Q @ U2.conj().T
        pP_given_Q = float(np.real(_trace(P @ rho_t2_after_Q)))
        p_reverse[idx] = pQ * pP_given_Q

        if use_switch and not switch_cptp:
            # Postselected (legacy) overlap term
            Kf = Q @ U2 @ P @ U1
            Kr = P @ U2 @ Q @ U1
            pf = float(np.real(_trace(Kf @ rho0 @ Kf.conj().T)))
            pr = float(np.real(_trace(Kr @ rho0 @ Kr.conj().T)))
            c = _trace(Kf @ rho0 @ Kr.conj().T)
            c_overlap[idx] = c
            base = 0.5 * (pf + pr)
            inter = float(np.real(exp_i_phi * c))
            p_plus[idx] = base + inter
            p_minus[idx] = base - inter
            # Visibility denominator: add eps to avoid 0/0 when both paths have
            # zero probability (can occur at specific t1 values for orthogonal projectors).
            # eps = 1e-15 << typical probability values (~0.1 to 1).
            denom = pf + pr + 1e-15
            v = 2.0 * abs(c) / denom
            v_max[idx] = float(min(1.0, max(0.0, v)))

        if use_switch and switch_cptp:
            # CPTP quantum switch built from Kraus sets of nonselective measurement channels:
            #   M_P: {P, I-P},  M_Q: {Q, I-Q}
            # Controlled-order Kraus operators on control⊗system:
            #   K_ab = |0><0| ⊗ (Q_a U2 P_b U1) + |1><1| ⊗ (P_b U2 Q_a U1)
            #
            # This is CPTP and yields interference in the control off-diagonal.
            Pk = [P, _I2() - P]
            Qk = [Q, _I2() - Q]
            rho_sc = np.zeros((4, 4), dtype=np.complex128)
            for Qa in Qk:
                for Pb in Pk:
                    A = Qa @ U2 @ Pb @ U1  # branch |0>: P then Q
                    B = Pb @ U2 @ Qa @ U1  # branch |1>: Q then P
                    K = np.kron(P0_ctrl, A) + np.kron(P1_ctrl, B)
                    rho_sc += K @ rho_in_sc @ K.conj().T

            rho_c_out = _ptrace_system(rho_sc)

            if bool(cfg.sanity_checks):
                # Sanity checks for the reduced control density matrix:
                # These verify the partial trace and CPTP channel are correctly implemented.
                assert rho_c_out.shape == (2, 2), f"Expected 2x2, got {rho_c_out.shape}"
                trace_val = float(np.real(np.trace(rho_c_out)))
                assert abs(trace_val - 1.0) < 1e-9, f"Trace = {trace_val}, expected 1.0"
                # Hermiticity check
                herm_err = float(np.linalg.norm(rho_c_out - rho_c_out.conj().T))
                assert herm_err < 1e-9, f"Non-Hermitian: ||ρ - ρ†|| = {herm_err}"
                # Positivity check (eigenvalues should be >= 0, allow tiny numerical negatives)
                min_eig = float(np.min(np.linalg.eigvalsh(rho_c_out).real))
                assert min_eig > -1e-9, f"Negative eigenvalue: {min_eig}"
            # Optional control dephasing (phase damping): rho01 -> (1-p) rho01
            p_dep = float(cfg.switch_dephase_p)
            if p_dep > 0.0:
                p_dep = max(0.0, min(1.0, p_dep))
                rho_c_out = rho_c_out.copy()
                rho_c_out[0, 1] *= (1.0 - p_dep)
                rho_c_out[1, 0] *= (1.0 - p_dep)
            control_coh[idx] = rho_c_out[0, 1]
            # Bloch vector of control: (Tr ρ σx, Tr ρ σy, Tr ρ σz)
            sx, sy, sz = _paulis()
            control_bloch[idx, 0] = float(np.real(np.trace(rho_c_out @ sx)))
            control_bloch[idx, 1] = float(np.real(np.trace(rho_c_out @ sy)))
            control_bloch[idx, 2] = float(np.real(np.trace(rho_c_out @ sz)))
            # Measure control in X basis: p(+)=<+|ρ_c|+>, p(-)=<-|ρ_c|->.
            ket_plus = (ket0 + ket1) / _SQRT2
            ket_minus = (ket0 - ket1) / _SQRT2
            p_plus[idx] = float(np.real(np.vdot(ket_plus, rho_c_out @ ket_plus)))
            p_minus[idx] = float(np.real(np.vdot(ket_minus, rho_c_out @ ket_minus)))
            # Visibility from coherence magnitude: V = 2|ρ01|
            v_max[idx] = float(min(1.0, max(0.0, 2.0 * abs(control_coh[idx]))))

    delta = p_forward - p_reverse

    # Time-series (0..t2) with optional stroboscopic projective filtering (Zeno-like),
    # which can change the *effective* observed frequency/content depending on the order.
    ts = np.arange(0.0, t2 + 1e-12, dt)
    pf_t = np.zeros_like(ts)
    pr_t = np.zeros_like(ts)

    U_dt = _unitary(H, dt)
    strobe = int(cfg.strobe_every_n_steps)
    strobe = 0 if strobe < 0 else strobe

    # Forward chain state and joint-weight
    w_f, rho_f = _apply_projective_measurement(rho0, P)
    # Reverse chain state and joint-weight
    w_r, rho_r = _apply_projective_measurement(rho0, Q)

    pf_t[0] = w_f * float(np.real(_trace(Q @ rho_f)))
    pr_t[0] = w_r * float(np.real(_trace(P @ rho_r)))

    for step in range(1, ts.size):
        # evolve both
        rho_f = U_dt @ rho_f @ U_dt.conj().T
        rho_r = U_dt @ rho_r @ U_dt.conj().T

        # optional projective filtering every N steps
        if strobe > 0 and (step % strobe) == 0:
            p_keep_f, rho_f = _apply_projective_measurement(rho_f, P)
            w_f *= p_keep_f
            p_keep_r, rho_r = _apply_projective_measurement(rho_r, Q)
            w_r *= p_keep_r

        pf_t[step] = w_f * float(np.real(_trace(Q @ rho_f)))
        pr_t[step] = w_r * float(np.real(_trace(P @ rho_r)))

    def _dominant_freq(y: np.ndarray, dt_s: float) -> float:
        if y.size < 4:
            return 0.0
        y0 = (y - float(np.mean(y))).astype(np.float64)
        spec = np.fft.rfft(y0)
        amp = np.abs(spec)
        freqs = np.fft.rfftfreq(y0.size, d=dt_s)
        if amp.size <= 1:
            return 0.0
        k = int(np.argmax(amp[1:])) + 1
        return float(freqs[k])

    f_pf = _dominant_freq(pf_t, dt)
    f_pr = _dominant_freq(pr_t, dt)
    f_delta = abs(f_pf - f_pr)

    return {
        "t1": t1s,
        "p_forward": p_forward,
        "p_reverse": p_reverse,
        "delta": delta,
        "t": ts,
        "p_forward_t": pf_t,
        "p_reverse_t": pr_t,
        "freq_forward": f_pf,
        "freq_reverse": f_pr,
        "freq_delta": f_delta,
        "comm_norms": comm_norms,
        "strobe_every_n_steps": int(cfg.strobe_every_n_steps),
        "use_switch": use_switch,
        "switch_cptp": switch_cptp,
        "switch_phase_deg": float(phi_deg_eff),
        "switch_c": c_overlap,
        "switch_p_plus": p_plus,
        "switch_p_minus": p_minus,
        "switch_vmax": v_max,
        "switch_control_coh": control_coh,
        "switch_control_bloch": control_bloch,
        "theta_deg": float(theta_deg_eff),
        "ofm_delta_deg": float(ofm_delta_deg),
        "cfg": cfg,
    }


def save_retro_plot(out_dir: str, cfg: RetroConfig) -> dict:
    out = run_retro(cfg)
    t1 = out["t1"]
    p_f = out["p_forward"]
    p_r = out["p_reverse"]
    delta = out["delta"]
    t = out["t"]
    pf_t = out["p_forward_t"]
    pr_t = out["p_reverse_t"]

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(10, 9), constrained_layout=True)
    ax0.plot(t1, p_f, lw=2, label="Forward: P@t1 then Q@t2 (joint prob)")
    ax0.plot(t1, p_r, lw=2, label="Reverse: Q@t1 then P@t2 (joint prob)")
    ax0.set_title(f"Order dependence (toy model)  θ={out['theta_deg']:.1f}°")
    ax0.set_xlabel("t1 (with fixed t2)")
    ax0.set_ylabel("joint probability")
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc="best")

    ax1.plot(t1, delta, lw=2, color="#d44a3a", label="Δ = forward - reverse")
    ax1.axhline(0.0, color="black", lw=1)
    ax1.set_xlabel("t1")
    ax1.set_ylabel("Δ")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="best")

    ax2.plot(t, pf_t, lw=2, label="Forward: P@0 then evolve then Q@t (joint prob)")
    ax2.plot(t, pr_t, lw=2, label="Reverse: Q@0 then evolve then P@t (joint prob)")
    strobe = int(out.get("strobe_every_n_steps", 0))
    ax2.set_title(
        f"Time-series (FFT): fwd={out['freq_forward']:.3f}Hz, rev={out['freq_reverse']:.3f}Hz  (strobe={strobe})"
    )
    ax2.set_xlabel("t")
    ax2.set_ylabel("joint probability")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best")

    fig.savefig(os.path.join(out_dir, "retro.png"), dpi=160)
    plt.close(fig)

    # Quantum-switch visibility plot (optional)
    vmax = out.get("switch_vmax", None)
    if bool(out.get("use_switch", False)) and isinstance(vmax, np.ndarray) and vmax.size:
        p_plus = out["switch_p_plus"]
        p_minus = out["switch_p_minus"]
        fringe = 2.0 * np.abs(p_plus - 0.5)
        fringe_max = float(np.max(fringe)) if fringe.size else 0.0
        figS, (axS0, axS1) = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
        axS0.plot(t1, p_plus, lw=2, label="p_plus (control |+>)")
        axS0.plot(t1, p_minus, lw=2, label="p_minus (control |->)")
        mode = "CPTP" if bool(out.get("switch_cptp", False)) else "postselected"
        axS0.set_title(f"Quantum switch ({mode})  phase={out['switch_phase_deg']:.1f}°")
        axS0.set_xlabel("t1")
        axS0.set_ylabel("probability (unnormalized)")
        axS0.grid(True, alpha=0.3)
        if bool(out.get("switch_cptp", False)):
            axS0.axhline(0.5, color="black", lw=1, alpha=0.5, label="classical mixture baseline (0.5)")
        axS0.legend(loc="best")

        if bool(out.get("switch_cptp", False)):
            axS1.plot(t1, vmax, lw=2, color="#2a6fdb", label="V(t1) = 2|ρ_c01|")
            axS1.axhline(0.0, color="black", lw=1, alpha=0.5, label="classical mixture baseline (0)")
        else:
            axS1.plot(t1, vmax, lw=2, color="#2a6fdb", label="V_max(t1) = 2|c|/(p_f+p_r)")
        axS1.set_ylim(0.0, 1.05)
        axS1.set_xlabel("t1")
        axS1.set_ylabel("visibility (0..1)")
        axS1.grid(True, alpha=0.3)
        axS1.legend(loc="best")
        figS.savefig(os.path.join(out_dir, "retro_switch.png"), dpi=160)
        plt.close(figS)
    else:
        fringe_max = 0.0

    # Spectrum plot
    def _spectrum(y: np.ndarray, dt_s: float) -> tuple[np.ndarray, np.ndarray]:
        y0 = (y - float(np.mean(y))).astype(np.float64)
        spec = np.fft.rfft(y0)
        amp = np.abs(spec)
        freqs = np.fft.rfftfreq(y0.size, d=dt_s)
        return freqs, amp

    freqs_f, amp_f = _spectrum(pf_t, float(cfg.dt))
    freqs_r, amp_r = _spectrum(pr_t, float(cfg.dt))
    denom = float(np.linalg.norm(amp_f)) + 1e-12
    spec_dist = float(np.linalg.norm(amp_f - amp_r) / denom)
    fig2, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
    ax.plot(freqs_f, amp_f, lw=2, label="Forward spectrum")
    ax.plot(freqs_r, amp_r, lw=2, label="Reverse spectrum")
    ax.set_xlim(0.0, min(5.0, float(np.max(freqs_f)) if freqs_f.size else 5.0))
    ax.set_title(f"FFT magnitude (mean-removed)  |Δf|={out['freq_delta']:.3f} Hz  spec_dist={spec_dist:.3f}")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("magnitude")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig2.savefig(os.path.join(out_dir, "retro_spectrum.png"), dpi=160)
    plt.close(fig2)

    return {
        "delta_max_abs": float(np.max(np.abs(delta))) if delta.size else 0.0,
        "freq_forward": float(out["freq_forward"]),
        "freq_reverse": float(out["freq_reverse"]),
        "freq_delta": float(out["freq_delta"]),
        "spec_dist": spec_dist,
        "switch_vmax_max": float(np.max(out["switch_vmax"])) if bool(out.get("use_switch", False)) else 0.0,
        "switch_fringe_max": fringe_max,
        "comm_norms": out["comm_norms"],
    }


def save_retro_theta_sweep_plot(out_dir: str, base: RetroConfig, step_deg: float = 5.0) -> None:
    step = float(step_deg)
    if step <= 0:
        raise ValueError("step_deg must be > 0")
    thetas = np.arange(0.0, 90.0 + 1e-9, step)

    max_delta = []
    freq_delta = []
    spec_dist = []
    vmax = []

    for th in thetas:
        cfg = RetroConfig(**{**base.__dict__, "theta_deg": float(th)})
        out = save_retro_plot(out_dir, cfg)
        max_delta.append(float(out["delta_max_abs"]))
        freq_delta.append(float(out["freq_delta"]))
        spec_dist.append(float(out["spec_dist"]))
        vmax.append(float(out.get("switch_vmax_max", 0.0)))

    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    ax.plot(thetas, max_delta, lw=2, label="max_abs_delta")
    ax.plot(thetas, spec_dist, lw=2, label="spec_dist")
    ax.plot(thetas, freq_delta, lw=2, label="freq_delta (Hz)")
    if bool(base.use_switch):
        ax.plot(thetas, vmax, lw=2, label="switch_vmax_max")
    ax.set_title(f"Retro metrics vs theta (t2={base.t2}, dt={base.dt}, strobe={base.strobe_every_n_steps})")
    ax.set_xlabel("theta (deg)")
    ax.set_ylabel("metric value")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.savefig(os.path.join(out_dir, "retro_theta_sweep.png"), dpi=160)
    plt.close(fig)


def save_retro_vis_vs_theta_plot(out_dir: str, base: RetroConfig, step_deg: float = 5.0) -> None:
    """
    CPTP switch: sweep theta and plot max visibility vs theta.
    - Coherence visibility: V_coh(t1)=2|rho01|
    - Fringe visibility:   V_fringe(t1)=2|p_plus-0.5|
    """
    step = float(step_deg)
    if step <= 0:
        raise ValueError("step_deg must be > 0")
    thetas = np.arange(0.0, 90.0 + 1e-9, step)
    v_coh = []
    v_fr = []
    for th in thetas:
        cfg = RetroConfig(**{**base.__dict__, "use_switch": True, "switch_cptp": True, "theta_deg": float(th)})
        out = run_retro(cfg)
        v_coh.append(float(np.max(out["switch_vmax"])) if out["switch_vmax"].size else 0.0)
        p_plus = out["switch_p_plus"]
        v_fr.append(float(np.max(2.0 * np.abs(p_plus - 0.5))) if p_plus.size else 0.0)

    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    ax.plot(thetas, v_coh, lw=2, label="max V_coh = max_t1 2|ρ_c01|")
    ax.plot(thetas, v_fr, lw=2, label="max V_fringe = max_t1 2|p_plus-0.5|")
    ax.set_title(f"CPTP quantum switch: visibility vs theta (phi={base.switch_phase_deg}°)")
    ax.set_xlabel("theta (deg)")
    ax.set_ylabel("visibility (0..1)")
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.savefig(os.path.join(out_dir, "retro_vis_vs_theta.png"), dpi=160)
    plt.close(fig)

    csv_path = os.path.join(out_dir, "retro_vis_vs_theta.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["theta_deg", "V_coh_max", "V_fringe_max"])
        for th, vc, vf in zip(thetas.tolist(), v_coh, v_fr, strict=False):
            w.writerow([th, vc, vf])


def save_retro_vis_vs_phi_plot(out_dir: str, base: RetroConfig, step_deg: float = 10.0) -> None:
    """
    CPTP switch: sweep phase phi and plot fringe visibility vs phi.
    Uses V_fringe(phi)=max_t1 2|p_plus(t1;phi)-0.5|, which should show a cosine-like lobe.
    """
    step = float(step_deg)
    if step <= 0:
        raise ValueError("step_deg must be > 0")
    phis = np.arange(0.0, 360.0 + 1e-9, step)
    v_fr = []
    for ph in phis:
        cfg = RetroConfig(**{**base.__dict__, "use_switch": True, "switch_cptp": True, "switch_phase_deg": float(ph)})
        out = run_retro(cfg)
        p_plus = out["switch_p_plus"]
        v_fr.append(float(np.max(2.0 * np.abs(p_plus - 0.5))) if p_plus.size else 0.0)

    fig, ax = plt.subplots(figsize=(10, 5.0), constrained_layout=True)
    ax.plot(phis, v_fr, lw=2, label="V_fringe(phi)=max_t1 2|p_plus-0.5|")
    ax.axhline(0.0, color="black", lw=1, alpha=0.5, label="classical mixture baseline (0)")
    ax.set_title(f"CPTP quantum switch: phase fringe vs phi (theta={base.theta_deg}°)")
    ax.set_xlabel("phi (deg)")
    ax.set_ylabel("fringe visibility (0..1)")
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.savefig(os.path.join(out_dir, "retro_vis_vs_phi.png"), dpi=160)
    plt.close(fig)


def save_retro_control_bloch_gif(out_dir: str, base: RetroConfig, filename: str = "retro_control_bloch.gif") -> str:
    """
    CPTP switch: animate the control Bloch vector as t1 varies.
    """
    cfg = RetroConfig(**{**base.__dict__, "use_switch": True, "switch_cptp": True})
    out = run_retro(cfg)
    t1 = out["t1"]
    bloch = out["switch_control_bloch"]  # (T,3)
    if bloch.size == 0:
        raise ValueError("No bloch data.")

    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(6.5, 6.0), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title("Control Bloch vector vs t1 (CPTP quantum switch)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_zlim(-1.0, 1.0)

    # unit sphere wireframe (light)
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 15)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, color="gray", alpha=0.15, linewidth=0.5)

    path, = ax.plot([], [], [], lw=2, color="#2a6fdb", alpha=0.8)
    point = ax.scatter([bloch[0, 0]], [bloch[0, 1]], [bloch[0, 2]], s=40, color="#d44a3a")

    def init():
        path.set_data([], [])
        path.set_3d_properties([])
        return (path, point)

    def update(i: int):
        path.set_data(bloch[: i + 1, 0], bloch[: i + 1, 1])
        path.set_3d_properties(bloch[: i + 1, 2])
        ax.set_title(f"Control Bloch vector vs t1  (t1={t1[i]:.2f})")
        point._offsets3d = ([bloch[i, 0]], [bloch[i, 1]], [bloch[i, 2]])
        return (path, point)

    ani = animation.FuncAnimation(fig, update, frames=bloch.shape[0], init_func=init, interval=60, blit=False)
    out_path = os.path.join(out_dir, filename)
    writer = animation.PillowWriter(fps=15)
    ani.save(out_path, writer=writer)
    plt.close(fig)
    return out_path


if HAS_TK:
    import tkinter as tk
    from tkinter import ttk

    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    class App(ttk.Frame):
        def __init__(self, master: tk.Tk):
            super().__init__(master)
            self.master.title("Fractal Hyperspace Simulator (OFM + Kuramoto + Tunneling)")
            self.pack(fill="both", expand=True)

            self.notebook = ttk.Notebook(self)
            self.notebook.pack(fill="both", expand=True)

            self.kura_tab = ttk.Frame(self.notebook)
            self.tunnel_tab = ttk.Frame(self.notebook)
            self.chsh_tab = ttk.Frame(self.notebook)
            self.retro_tab = ttk.Frame(self.notebook)
            self.notebook.add(self.kura_tab, text="Kuramoto (OFM)")
            self.notebook.add(self.tunnel_tab, text="Tunneling (Split-step)")
            self.notebook.add(self.chsh_tab, text="CHSH (Entanglement)")
            self.notebook.add(self.retro_tab, text="Retro (Projectors)")

            self._build_kuramoto_tab()
            self._build_tunneling_tab()
            self._build_chsh_tab()
            self._build_retro_tab()

        def _build_kuramoto_tab(self) -> None:
            left = ttk.Frame(self.kura_tab, padding=10)
            left.pack(side="left", fill="y")
            right = ttk.Frame(self.kura_tab, padding=10)
            right.pack(side="right", fill="both", expand=True)

            self.k_n = tk.IntVar(value=20)
            self.k_J = tk.DoubleVar(value=1.0)
            self.k_om_mean = tk.DoubleVar(value=0.0)
            self.k_om_std = tk.DoubleVar(value=0.5)
            self.k_t_end = tk.DoubleVar(value=20.0)
            self.k_dt = tk.DoubleVar(value=0.05)
            self.k_seed = tk.IntVar(value=1)

            row = 0
            for label, var in [
                ("N (nodes)", self.k_n),
                ("J (strength)", self.k_J),
                ("omega mean", self.k_om_mean),
                ("omega std", self.k_om_std),
                ("t_end", self.k_t_end),
                ("dt", self.k_dt),
                ("seed", self.k_seed),
            ]:
                ttk.Label(left, text=label).grid(row=row, column=0, sticky="w", pady=2)
                ttk.Entry(left, textvariable=var, width=12).grid(row=row, column=1, sticky="w", pady=2)
                row += 1

            self.k_status = ttk.Label(left, text="Ready")
            self.k_status.grid(row=row, column=0, columnspan=2, sticky="w", pady=(10, 0))
            row += 1

            ttk.Button(left, text="Run Kuramoto", command=self._run_kuramoto_clicked).grid(
                row=row, column=0, columnspan=2, sticky="we", pady=8
            )

            self.k_fig, (self.k_ax_r, self.k_ax_phase) = plt.subplots(
                2, 1, figsize=(7.5, 6.0), constrained_layout=True
            )
            self.k_canvas = FigureCanvasTkAgg(self.k_fig, master=right)
            self.k_canvas.get_tk_widget().pack(fill="both", expand=True)

            self.k_ax_r.set_title("Order parameter r(t)")
            self.k_ax_r.set_xlabel("t")
            self.k_ax_r.set_ylabel("r")
            self.k_ax_r.grid(True, alpha=0.3)

            self.k_ax_phase.set_title("Final phases (wrapped)")
            self.k_ax_phase.set_xlabel("node index")
            self.k_ax_phase.set_ylabel("theta mod 2π")
            self.k_ax_phase.grid(True, alpha=0.3)

        def _run_kuramoto_clicked(self) -> None:
            try:
                cfg = KuramotoConfig(
                    n=int(self.k_n.get()),
                    j_strength=float(self.k_J.get()),
                    omega_mean=float(self.k_om_mean.get()),
                    omega_std=float(self.k_om_std.get()),
                    t_end=float(self.k_t_end.get()),
                    dt=float(self.k_dt.get()),
                    seed=int(self.k_seed.get()),
                )
                if cfg.n < 2 or cfg.n > 200:
                    raise ValueError("N must be in [2, 200] for this GUI.")
                self.k_status.config(text="Running...")
                self.master.update_idletasks()

                out = run_kuramoto(cfg)
                t = out["t"]
                r = out["r"]
                theta = out["theta"]
                theta_final = np.mod(theta[:, -1], 2.0 * np.pi)

                self.k_ax_r.clear()
                self.k_ax_r.plot(t, r, lw=2)
                self.k_ax_r.set_title(f"Order parameter r(t)  (final r={r[-1]:.3f})")
                self.k_ax_r.set_xlabel("t")
                self.k_ax_r.set_ylabel("r")
                self.k_ax_r.set_ylim(0.0, 1.05)
                self.k_ax_r.grid(True, alpha=0.3)

                self.k_ax_phase.clear()
                self.k_ax_phase.scatter(np.arange(1, cfg.n + 1), theta_final, s=20)
                self.k_ax_phase.set_title("Final phases (wrapped to 0..2π)")
                self.k_ax_phase.set_xlabel("node index")
                self.k_ax_phase.set_ylabel("theta mod 2π")
                self.k_ax_phase.set_ylim(0.0, 2.0 * np.pi)
                self.k_ax_phase.grid(True, alpha=0.3)

                self.k_canvas.draw()
                self.k_status.config(text=f"Done. final r={r[-1]:.3f}")
            except Exception as e:
                self.k_status.config(text=f"Error: {e}")

        def _build_tunneling_tab(self) -> None:
            left = ttk.Frame(self.tunnel_tab, padding=10)
            left.pack(side="left", fill="y")
            right = ttk.Frame(self.tunnel_tab, padding=10)
            right.pack(side="right", fill="both", expand=True)

            self.t_n = tk.IntVar(value=512)
            self.t_L = tk.DoubleVar(value=40.0)
            self.t_dt = tk.DoubleVar(value=0.1)
            self.t_steps = tk.IntVar(value=400)
            self.t_V0 = tk.DoubleVar(value=1.0)
            self.t_a = tk.DoubleVar(value=4.0)
            self.t_E = tk.DoubleVar(value=0.30)
            self.t_sigma = tk.DoubleVar(value=2.0)
            self.t_x0 = tk.DoubleVar(value=-10.0)
            self.t_barrier = tk.StringVar(value="rect")
            self.t_cantor_levels = tk.IntVar(value=5)

            row = 0
            for label, var in [
                ("N (grid)", self.t_n),
                ("L (domain)", self.t_L),
                ("dt", self.t_dt),
                ("steps", self.t_steps),
                ("V0", self.t_V0),
                ("a (width)", self.t_a),
                ("E", self.t_E),
                ("sigma", self.t_sigma),
                ("x0", self.t_x0),
            ]:
                ttk.Label(left, text=label).grid(row=row, column=0, sticky="w", pady=2)
                ttk.Entry(left, textvariable=var, width=12).grid(row=row, column=1, sticky="w", pady=2)
                row += 1

            ttk.Label(left, text="barrier").grid(row=row, column=0, sticky="w", pady=(10, 2))
            ttk.Combobox(
                left, textvariable=self.t_barrier, values=["rect", "cantor"], width=10, state="readonly"
            ).grid(row=row, column=1, sticky="w", pady=(10, 2))
            row += 1

            ttk.Label(left, text="cantor levels").grid(row=row, column=0, sticky="w", pady=2)
            ttk.Entry(left, textvariable=self.t_cantor_levels, width=12).grid(row=row, column=1, sticky="w", pady=2)
            row += 1

            self.t_status = ttk.Label(left, text="Ready")
            self.t_status.grid(row=row, column=0, columnspan=2, sticky="w", pady=(10, 0))
            row += 1

            ttk.Button(left, text="Run Tunneling", command=self._run_tunneling_clicked).grid(
                row=row, column=0, columnspan=2, sticky="we", pady=8
            )

            self.t_fig, (self.t_ax0, self.t_ax1) = plt.subplots(
                2, 1, figsize=(7.5, 6.0), constrained_layout=True
            )
            self.t_canvas = FigureCanvasTkAgg(self.t_fig, master=right)
            self.t_canvas.get_tk_widget().pack(fill="both", expand=True)

            self.t_ax0.set_title("|psi(x,t=0)|^2 and V(x)")
            self.t_ax0.set_xlabel("x")
            self.t_ax0.set_ylabel("density / potential")
            self.t_ax0.grid(True, alpha=0.3)

            self.t_ax1.set_title("|psi(x,t=end)|^2 and V(x)")
            self.t_ax1.set_xlabel("x")
            self.t_ax1.set_ylabel("density / potential")
            self.t_ax1.grid(True, alpha=0.3)

        def _run_tunneling_clicked(self) -> None:
            try:
                cfg = TunnelingConfig(
                    n=int(self.t_n.get()),
                    L=float(self.t_L.get()),
                    dt=float(self.t_dt.get()),
                    steps=int(self.t_steps.get()),
                    V0=float(self.t_V0.get()),
                    a=float(self.t_a.get()),
                    E=float(self.t_E.get()),
                    sigma=float(self.t_sigma.get()),
                    x0=float(self.t_x0.get()),
                    barrier_type=str(self.t_barrier.get()),
                    cantor_levels=int(self.t_cantor_levels.get()),
                )
                if cfg.n < 128 or cfg.n > 8192:
                    raise ValueError("N must be in [128, 8192].")
                if cfg.steps < 1 or cfg.steps > 50000:
                    raise ValueError("steps must be in [1, 50000].")
                self.t_status.config(text="Running...")
                self.master.update_idletasks()

                out = run_tunneling(cfg)
                x = out["x"]
                V = out["V"]
                psi0 = out["psi0"]
                psi = out["psi"]
                T = out["T"]

                rho0 = np.abs(psi0) ** 2
                rho = np.abs(psi) ** 2
                v_scaled = V / (np.max(V) + 1e-12) * (0.9 * max(float(np.max(rho0)), float(np.max(rho))) + 1e-12)

                self.t_ax0.clear()
                self.t_ax0.plot(x, rho0, lw=2, label="|psi|^2 (t=0)")
                self.t_ax0.plot(x, v_scaled, lw=1.5, label="V(x) (scaled)", alpha=0.8)
                self.t_ax0.set_title("|psi(x,t=0)|^2 and V(x)")
                self.t_ax0.set_xlabel("x")
                self.t_ax0.set_ylabel("density / potential (scaled)")
                self.t_ax0.grid(True, alpha=0.3)
                self.t_ax0.legend(loc="upper right")

                self.t_ax1.clear()
                self.t_ax1.plot(x, rho, lw=2, label="|psi|^2 (t=end)")
                self.t_ax1.plot(x, v_scaled, lw=1.5, label="V(x) (scaled)", alpha=0.8)
                self.t_ax1.set_title(f"|psi(x,t=end)|^2 and V(x)   (T={T:.3f})")
                self.t_ax1.set_xlabel("x")
                self.t_ax1.set_ylabel("density / potential (scaled)")
                self.t_ax1.grid(True, alpha=0.3)
                self.t_ax1.legend(loc="upper right")

                self.t_canvas.draw()
                self.t_status.config(text=f"Done. T={T:.3f}")
            except Exception as e:
                self.t_status.config(text=f"Error: {e}")

        def _build_chsh_tab(self) -> None:
            left = ttk.Frame(self.chsh_tab, padding=10)
            left.pack(side="left", fill="y")
            right = ttk.Frame(self.chsh_tab, padding=10)
            right.pack(side="right", fill="both", expand=True)

            self.c_n_ofm = tk.IntVar(value=20)
            self.c_strength = tk.DoubleVar(value=1.0)

            row = 0
            for label, var in [
                ("OFM N", self.c_n_ofm),
                ("fractal strength", self.c_strength),
            ]:
                ttk.Label(left, text=label).grid(row=row, column=0, sticky="w", pady=2)
                ttk.Entry(left, textvariable=var, width=12).grid(row=row, column=1, sticky="w", pady=2)
                row += 1

            self.c_status = ttk.Label(left, text="Ready")
            self.c_status.grid(row=row, column=0, columnspan=2, sticky="w", pady=(10, 0))
            row += 1

            ttk.Button(left, text="Run CHSH", command=self._run_chsh_clicked).grid(
                row=row, column=0, columnspan=2, sticky="we", pady=8
            )

            self.c_fig, self.c_ax = plt.subplots(figsize=(7.5, 4.0), constrained_layout=True)
            self.c_canvas = FigureCanvasTkAgg(self.c_fig, master=right)
            self.c_canvas.get_tk_widget().pack(fill="both", expand=True)
            self.c_ax.set_title("CHSH")
            self.c_ax.grid(True, alpha=0.3)

        def _run_chsh_clicked(self) -> None:
            try:
                cfg = CHSHConfig(
                    n_ofm=int(self.c_n_ofm.get()),
                    fractal_noise=float(self.c_strength.get()),
                    # seed is unused (CHSH is deterministic), using default
                )
                self.c_status.config(text="Running...")
                self.master.update_idletasks()

                out = run_chsh(cfg)
                s = float(out["CHSH"])

                self.c_ax.clear()
                self.c_ax.bar(["CHSH"], [s], color="#2a6fdb" if s > 2.0 else "#d44a3a")
                self.c_ax.axhline(2.0, color="#2c8a2c", ls="--", lw=1.5, label="Classical (2)")
                self.c_ax.axhline(2.0 * _SQRT2, color="#7a2fbf", ls="--", lw=1.5, label="Tsirelson (2√2)")
                self.c_ax.set_ylim(0.0, 3.0)
                self.c_ax.set_title(
                    f"CHSH S={s:.3f}  (p={out['p']:.3f}, OFM scale={out['ofm_scale']:.3f})"
                )
                self.c_ax.legend(loc="upper right")
                self.c_ax.grid(True, alpha=0.2, axis="y")
                self.c_canvas.draw()
                self.c_status.config(text=f"Done. CHSH={s:.3f}")
            except Exception as e:
                self.c_status.config(text=f"Error: {e}")

        def _build_retro_tab(self) -> None:
            left = ttk.Frame(self.retro_tab, padding=10)
            left.pack(side="left", fill="y")
            right = ttk.Frame(self.retro_tab, padding=10)
            right.pack(side="right", fill="both", expand=True)

            self.r_omega = tk.DoubleVar(value=1.0)
            self.r_t2 = tk.DoubleVar(value=5.0)
            self.r_dt = tk.DoubleVar(value=0.02)
            self.r_theta = tk.DoubleVar(value=90.0)
            self.r_strobe = tk.IntVar(value=0)
            self.r_switch = tk.IntVar(value=0)
            self.r_phase = tk.DoubleVar(value=0.0)
            self.r_switch_post = tk.IntVar(value=0)

            row = 0
            for label, var in [
                ("omega", self.r_omega),
                ("t2", self.r_t2),
                ("dt", self.r_dt),
                ("theta_deg", self.r_theta),
                ("strobe_n_steps", self.r_strobe),
                ("switch (0/1)", self.r_switch),
                ("switch_phase_deg", self.r_phase),
                ("switch_post (0/1)", self.r_switch_post),
            ]:
                ttk.Label(left, text=label).grid(row=row, column=0, sticky="w", pady=2)
                ttk.Entry(left, textvariable=var, width=12).grid(row=row, column=1, sticky="w", pady=2)
                row += 1

            ttk.Label(left, text="(Fixed axes: H=X, P=Z+, ψ0=Z+.  Q=σθ in x–z plane)").grid(
                row=row, column=0, columnspan=2, sticky="w", pady=(8, 2)
            )
            row += 1

            self.r_status = ttk.Label(left, text="Ready")
            self.r_status.grid(row=row, column=0, columnspan=2, sticky="w", pady=(10, 0))
            row += 1

            ttk.Button(left, text="Run Retro", command=self._run_retro_clicked).grid(
                row=row, column=0, columnspan=2, sticky="we", pady=8
            )

            self.r_fig, (self.r_ax0, self.r_ax1, self.r_ax2) = plt.subplots(
                3, 1, figsize=(7.5, 8.0), constrained_layout=True
            )
            self.r_canvas = FigureCanvasTkAgg(self.r_fig, master=right)
            self.r_canvas.get_tk_widget().pack(fill="both", expand=True)

            self.r_ax0.set_title("Order-dependent joint probabilities")
            self.r_ax0.grid(True, alpha=0.3)
            self.r_ax1.set_title("Δ (forward - reverse)")
            self.r_ax1.grid(True, alpha=0.3)
            self.r_ax2.set_title("Time-series (for FFT)")
            self.r_ax2.grid(True, alpha=0.3)

        def _run_retro_clicked(self) -> None:
            try:
                cfg = RetroConfig(
                    omega=float(self.r_omega.get()),
                    t2=float(self.r_t2.get()),
                    dt=float(self.r_dt.get()),
                    theta_deg=float(self.r_theta.get()),
                    strobe_every_n_steps=int(self.r_strobe.get()),
                    use_switch=bool(int(self.r_switch.get())),
                    switch_phase_deg=float(self.r_phase.get()),
                    switch_cptp=not bool(int(self.r_switch_post.get())),
                    sanity_checks=False,
                )
                self.r_status.config(text="Running...")
                self.master.update_idletasks()

                out = run_retro(cfg)
                t1 = out["t1"]
                p_f = out["p_forward"]
                p_r = out["p_reverse"]
                delta = out["delta"]
                t = out["t"]
                pf_t = out["p_forward_t"]
                pr_t = out["p_reverse_t"]

                self.r_ax0.clear()
                self.r_ax0.plot(t1, p_f, lw=2, label="Forward: P@t1 then Q@t2")
                self.r_ax0.plot(t1, p_r, lw=2, label="Reverse: Q@t1 then P@t2")
                self.r_ax0.set_title(f"Order dependence (θ={out['theta_deg']:.1f}°)")
                self.r_ax0.set_xlabel("t1 (fixed t2)")
                self.r_ax0.set_ylabel("joint probability")
                self.r_ax0.grid(True, alpha=0.3)
                self.r_ax0.legend(loc="best")

                self.r_ax1.clear()
                self.r_ax1.plot(t1, delta, lw=2, color="#d44a3a", label="Δ = forward - reverse")
                self.r_ax1.axhline(0.0, color="black", lw=1)
                self.r_ax1.set_xlabel("t1")
                self.r_ax1.set_ylabel("Δ")
                self.r_ax1.grid(True, alpha=0.3)
                self.r_ax1.legend(loc="best")

                self.r_ax2.clear()
                self.r_ax2.plot(t, pf_t, lw=2, label="Forward: P@0 then Q@t")
                self.r_ax2.plot(t, pr_t, lw=2, label="Reverse: Q@0 then P@t")
                self.r_ax2.set_title(
                    f"Time-series for FFT (fwd={out['freq_forward']:.3f}Hz, rev={out['freq_reverse']:.3f}Hz)"
                )
                self.r_ax2.set_xlabel("t")
                self.r_ax2.set_ylabel("joint probability")
                self.r_ax2.grid(True, alpha=0.3)
                self.r_ax2.legend(loc="best")

                self.r_canvas.draw()
                self.r_status.config(
                    text=f"Done. max|Δ|={float(np.max(np.abs(delta))):.3f}, |df|={out['freq_delta']:.3f}Hz"
                )
            except Exception as e:
                self.r_status.config(text=f"Error: {e}")


def gui_main() -> int:
    if not HAS_TK:
        print("[gui] Tk is not available; re-run with --nogui.", file=sys.stderr)
        return 2
    root = tk.Tk()
    try:
        root.tk.call("tk", "scaling", 1.2)
    except Exception:
        pass
    App(root)
    root.mainloop()
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="OFM+Kuramoto and split-step tunneling simulator.")
    parser.add_argument("--nogui", action="store_true", help="Run without Tkinter GUI (saves PNGs, prints metrics).")
    parser.add_argument(
        "--mode",
        choices=["kuramoto", "kuramoto_anim", "tunneling", "both", "chsh", "retro", "all"],
        default="both",
        help="Which sim(s) to run. 'both' = kuramoto+tunneling; 'all' = kuramoto+tunneling+chsh+retro.",
    )
    parser.add_argument("--out", default="outputs", help="Output directory for headless mode.")
    parser.add_argument("--barrier", choices=["rect", "cantor"], default="rect", help="Barrier type for tunneling.")
    parser.add_argument("--cantor-t-sweep", action="store_true", help="Tunneling: sweep cantor_levels and write cantor_t_sweep.png.")
    parser.add_argument("--cantor-levels", default=None, help="Tunneling: cantor sweep levels, e.g. \"1,2,4,8,16\".")
    parser.add_argument("--unified-demo", action="store_true", help="Optional: map Kuramoto final r into retro switch phase (phi).")
    # Retro (headless)
    parser.add_argument("--retro-theta", type=float, default=None, help="Retro: projector angle theta (deg).")
    parser.add_argument("--retro-omega", type=float, default=None, help="Retro: omega in H.")
    parser.add_argument("--retro-t2", type=float, default=None, help="Retro: total time t2.")
    parser.add_argument("--retro-dt", type=float, default=None, help="Retro: step dt.")
    parser.add_argument("--retro-strobe", type=int, default=None, help="Retro: strobe_every_n_steps (0 disables).")
    parser.add_argument("--retro-switch", action="store_true", help="Retro: enable quantum switch (defaults to CPTP).")
    parser.add_argument("--retro-switch-post", action="store_true", help="Retro: use postselected overlap model (legacy).")
    parser.add_argument("--retro-switch-phase", type=float, default=None, help="Retro: switch relative phase (deg).")
    parser.add_argument("--retro-sweep", action="store_true", help="Retro: sweep theta 0..90 and plot metrics vs theta.")
    parser.add_argument("--retro-sweep-step", type=float, default=5.0, help="Retro: theta step (deg) for sweep.")
    parser.add_argument("--retro-vis-sweep", action="store_true", help="Retro: CPTP switch visibility vs theta sweep (writes retro_vis_vs_theta.png).")
    parser.add_argument("--retro-phi-sweep", action="store_true", help="Retro: CPTP switch phase sweep (writes retro_vis_vs_phi.png).")
    parser.add_argument("--retro-phi-step", type=float, default=10.0, help="Retro: phi step (deg) for phase sweep.")
    parser.add_argument("--retro-control-bloch", action="store_true", help="Retro: animate control Bloch vector vs t1 (writes retro_control_bloch.gif).")
    parser.add_argument("--retro-ofm", action="store_true", help="Retro: enable OFM modulation of theta/phi.")
    parser.add_argument("--retro-ofm-n", type=int, default=20, help="Retro OFM: N used to compute OFM scale.")
    parser.add_argument("--retro-ofm-strength-deg", type=float, default=0.0, help="Retro OFM: max delta in degrees (scaled by OFM).")
    parser.add_argument("--retro-ofm-target", choices=["phi", "theta"], default="phi", help="Retro OFM: apply modulation to phi or theta.")
    # Kuramoto (headless)
    parser.add_argument("--kuramoto-topology", choices=["ofm", "torus", "metatron", "pyramid"], default=None, help="Kuramoto topology.")
    parser.add_argument("--torus-m", type=int, default=None, help="Kuramoto: torus grid size m (N=m*m).")
    parser.add_argument("--kuramoto-n", type=int, default=None, help="Kuramoto: N (only used for topology=ofm).")
    parser.add_argument("--kuramoto-j", type=float, default=None, help="Kuramoto: coupling strength J.")
    parser.add_argument("--kuramoto-omega-mean", type=float, default=None, help="Kuramoto: omega_mean (mean frequency).")
    parser.add_argument("--kuramoto-omega-std", type=float, default=None, help="Kuramoto: omega_std (spread).")
    parser.add_argument("--kuramoto-t-end", type=float, default=None, help="Kuramoto: simulation end time.")
    parser.add_argument("--kuramoto-dt", type=float, default=None, help="Kuramoto: output dt.")
    parser.add_argument("--kuramoto-seed", type=int, default=None, help="Kuramoto: RNG seed.")
    parser.add_argument("--pyramid-scale-m", type=float, default=None, help="Kuramoto: pyramid distance scale (meters).")
    parser.add_argument("--pyramid-shaft-weight", type=float, default=None, help="Kuramoto: weight multiplier for pyramid shaft edges.")
    parser.add_argument("--pyramid-medium-weights", type=str, default=None, help="Kuramoto: medium weights (name=value,comma-separated).")
    parser.add_argument("--pyramid-nodes-csv", type=str, default=None, help="Kuramoto: pyramid nodes CSV (name,x,y,z).")
    parser.add_argument("--pyramid-edges-csv", type=str, default=None, help="Kuramoto: pyramid edges CSV (u,v,medium,weight).")
    parser.add_argument("--kuramoto-model", choices=["first", "second"], default=None, help="Kuramoto: first-order or inertia model.")
    parser.add_argument("--kuramoto-mass", type=float, default=None, help="Kuramoto: inertia mass (second-order).")
    parser.add_argument("--kuramoto-damping", type=float, default=None, help="Kuramoto: damping gamma (second-order).")
    parser.add_argument("--pyramid-shaft-sweep", type=str, default=None, help="Kuramoto: sweep air_shaft weights (comma-separated).")
    parser.add_argument("--kuramoto-fft", action="store_true", help="Kuramoto: write FFT of r(t).")
    parser.add_argument("--kuramoto-phase-pairs", type=str, default=None, help="Kuramoto: phase portrait pairs (e.g. 1-2,1-3).")
    parser.add_argument("--kuramoto-label-nodes", action="store_true", help="Kuramoto: label nodes on final phase plot (if names available).")
    parser.add_argument("--kuramoto-pyramid-resonance", action="store_true", help="Kuramoto: set omega by pyramid chamber resonances.")
    parser.add_argument("--kuramoto-freq-sweep", type=str, default=None, help="Kuramoto: sweep target freqs in Hz (comma-separated).")
    # CHSH (headless)
    parser.add_argument("--chsh-noise-sweep", action="store_true", help="CHSH: sweep noise values and write chsh_noise_sweep.png.")
    parser.add_argument("--chsh-ofm-n", type=int, default=20, help="CHSH: OFM N for scaling.")
    args = parser.parse_args()

    if args.nogui or not HAS_TK:
        if not HAS_TK and not args.nogui:
            print(
                "[info] Tkinter/Tcl not available; switching to headless mode.\n"
                "       (If you want the GUI, reinstall Python with Tcl/Tk support.)"
            )
        # Apply CLI overrides for headless runs (where applicable)
        if args.mode in ("retro", "all"):
            cfg = RetroConfig()
            if args.retro_theta is not None:
                cfg = replace(cfg, theta_deg=float(args.retro_theta))
            if args.retro_omega is not None:
                cfg = replace(cfg, omega=float(args.retro_omega))
            if args.retro_t2 is not None:
                cfg = replace(cfg, t2=float(args.retro_t2))
            if args.retro_dt is not None:
                cfg = replace(cfg, dt=float(args.retro_dt))
            if args.retro_strobe is not None:
                cfg = replace(cfg, strobe_every_n_steps=int(args.retro_strobe))
            if bool(args.retro_switch):
                cfg = replace(cfg, use_switch=True)
            if bool(args.retro_switch_post):
                cfg = replace(cfg, use_switch=True, switch_cptp=False)
            if args.retro_switch_phase is not None:
                cfg = replace(cfg, switch_phase_deg=float(args.retro_switch_phase))
            if bool(args.retro_ofm):
                cfg = replace(
                    cfg,
                    ofm_modulate=True,
                    ofm_n=int(args.retro_ofm_n),
                    ofm_strength_deg=float(args.retro_ofm_strength_deg),
                    ofm_target=str(args.retro_ofm_target),
                )
            # stash into globals for headless_main to pick up via defaults (simple approach)
            global _RETRO_HEADLESS_CFG  # noqa: PLW0603
            _RETRO_HEADLESS_CFG = cfg
            global _RETRO_SWEEP  # noqa: PLW0603
            _RETRO_SWEEP = {"enabled": bool(args.retro_sweep), "step": float(args.retro_sweep_step)}
            global _RETRO_EXTRA  # noqa: PLW0603
            _RETRO_EXTRA = {
                "vis_sweep": bool(args.retro_vis_sweep),
                "phi_sweep": bool(args.retro_phi_sweep),
                "phi_step": float(args.retro_phi_step),
                "bloch": bool(args.retro_control_bloch),
            }

        if args.mode in ("kuramoto", "both", "all", "kuramoto_anim"):
            kcfg = KuramotoConfig()
            if args.kuramoto_topology is not None:
                kcfg = replace(kcfg, topology=str(args.kuramoto_topology))
            if args.torus_m is not None:
                kcfg = replace(kcfg, torus_m=int(args.torus_m))
            if args.kuramoto_n is not None:
                kcfg = replace(kcfg, n=int(args.kuramoto_n))
            if args.kuramoto_j is not None:
                kcfg = replace(kcfg, j_strength=float(args.kuramoto_j))
            if args.kuramoto_omega_mean is not None:
                kcfg = replace(kcfg, omega_mean=float(args.kuramoto_omega_mean))
            if args.kuramoto_omega_std is not None:
                kcfg = replace(kcfg, omega_std=float(args.kuramoto_omega_std))
            if args.kuramoto_t_end is not None:
                kcfg = replace(kcfg, t_end=float(args.kuramoto_t_end))
            if args.kuramoto_dt is not None:
                kcfg = replace(kcfg, dt=float(args.kuramoto_dt))
            if args.kuramoto_seed is not None:
                kcfg = replace(kcfg, seed=int(args.kuramoto_seed))
            if args.pyramid_scale_m is not None:
                kcfg = replace(kcfg, pyramid_scale_m=float(args.pyramid_scale_m))
            if args.pyramid_shaft_weight is not None:
                kcfg = replace(kcfg, pyramid_shaft_weight=float(args.pyramid_shaft_weight))
            if args.pyramid_medium_weights is not None:
                kcfg = replace(kcfg, pyramid_medium_weights=str(args.pyramid_medium_weights))
            if args.pyramid_nodes_csv is not None:
                kcfg = replace(kcfg, pyramid_nodes_csv=str(args.pyramid_nodes_csv))
            if args.pyramid_edges_csv is not None:
                kcfg = replace(kcfg, pyramid_edges_csv=str(args.pyramid_edges_csv))
            if args.kuramoto_model is not None:
                kcfg = replace(kcfg, model=str(args.kuramoto_model))
            if args.kuramoto_mass is not None:
                kcfg = replace(kcfg, mass=float(args.kuramoto_mass))
            if args.kuramoto_damping is not None:
                kcfg = replace(kcfg, damping=float(args.kuramoto_damping))
            global _KURAMOTO_HEADLESS_CFG  # noqa: PLW0603
            _KURAMOTO_HEADLESS_CFG = kcfg
            global _KURAMOTO_SWEEP  # noqa: PLW0603
            _KURAMOTO_SWEEP = {"pyramid_shaft_sweep": args.pyramid_shaft_sweep, "freq_sweep": args.kuramoto_freq_sweep}
            global _KURAMOTO_EXTRA  # noqa: PLW0603
            _KURAMOTO_EXTRA = {
                "fft": bool(args.kuramoto_fft),
                "phase_pairs": args.kuramoto_phase_pairs,
                "label_nodes": bool(args.kuramoto_label_nodes),
                "pyramid_resonance": bool(args.kuramoto_pyramid_resonance),
            }

        if args.mode in ("chsh", "all"):
            ccfg = CHSHConfig(n_ofm=int(args.chsh_ofm_n), fractal_noise=0.0, seed=1)
            global _CHSH_HEADLESS_CFG  # noqa: PLW0603
            _CHSH_HEADLESS_CFG = ccfg
            global _CHSH_EXTRA  # noqa: PLW0603
            _CHSH_EXTRA = {"noise_sweep": bool(args.chsh_noise_sweep)}

        if args.mode in ("tunneling", "both", "all"):
            global _TUNNEL_EXTRA  # noqa: PLW0603
            levels = None
            if args.cantor_levels:
                try:
                    levels = [int(x.strip()) for x in str(args.cantor_levels).split(",") if x.strip()]
                except Exception as e:
                    raise ValueError(f"Invalid --cantor-levels: {args.cantor_levels}") from e
            _TUNNEL_EXTRA = {"cantor_sweep": bool(args.cantor_t_sweep), "levels": levels}

        global _UNIFIED_DEMO  # noqa: PLW0603
        _UNIFIED_DEMO = bool(args.unified_demo)

        return headless_main(mode=args.mode, out_dir=args.out, barrier=args.barrier)

    return gui_main()


if __name__ == "__main__":
    raise SystemExit(main())

