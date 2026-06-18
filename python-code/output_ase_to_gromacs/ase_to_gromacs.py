##  #!/usr/bin/env python3

"""
Convert an ASE .traj file (from an NPT-MD simulation) to a GROMACS .trr file.

Dependencies:
    pip install ase MDAnalysis

Usage:
    python traj_to_trr.py input.traj output.trr [--dt 1.0] [--skip 1]
"""

import argparse
import numpy as np

try:
    from ase.io.trajectory import Trajectory
except ImportError:
    raise ImportError("ASE is required: pip install ase")

try:
    import MDAnalysis as mda
    from MDAnalysis.coordinates.TRR import TRRWriter
except ImportError:
    raise ImportError("MDAnalysis is required: pip install MDAnalysis")


def traj_to_trr(
    input_traj: str,
    output_trr: str,
    dt: float = 1.0,
    skip: int = 1,
) -> None:
    """
    Convert an ASE Trajectory file to a GROMACS TRR file.

    Parameters
    ----------
    input_traj : str
        Path to the ASE .traj file.
    output_trr : str
        Path for the output .trr file.
    dt : float
        Timestep between saved frames in picoseconds (default: 1.0 ps).
        Used to assign simulation time to each frame.
    skip : int
        Write every Nth frame (default: 1 = all frames).
    """
    print(f"Reading ASE trajectory: {input_traj}")
    traj = Trajectory(input_traj, "r")

    # Read all (or strided) frames into memory first so we know the count
    frames = [traj[i] for i in range(0, len(traj), skip)]
    n_frames = len(frames)
    n_atoms = len(frames[0])

    print(f"  Frames to write : {n_frames}  (skip={skip})")
    print(f"  Atoms per frame : {n_atoms}")
    print(f"  Timestep        : {dt} ps")

    # ------------------------------------------------------------------
    # Build a minimal MDAnalysis Universe so TRRWriter has atom metadata.
    # We only need positions/velocities/forces – element labels are enough.
    # ------------------------------------------------------------------
    symbols = frames[0].get_chemical_symbols()

    # ------------------------------------------------------------------
    # Probe the first frame to know which optional arrays are present.
    # MDAnalysis Timestep must be pre-allocated with the right flags or
    # assigning .velocities / .forces raises NoDataError.
    # ------------------------------------------------------------------
    _f0 = frames[0]
    _has_vel = _f0.get_velocities() is not None
    _has_frc = True
    try:
        _f0.get_forces(apply_constraint=False)
    except Exception:
        _has_frc = False

    print(f"  Has velocities  : {_has_vel}")
    print(f"  Has forces      : {_has_frc}")

    u = mda.Universe.empty(
        n_atoms,
        n_residues=1,
        atom_resindex=np.zeros(n_atoms, dtype=int),
        residue_segindex=np.zeros(1, dtype=int),
        trajectory=True,
    )
    u.add_TopologyAttr("name", symbols)
    u.add_TopologyAttr("resname", ["SYS"])
    u.add_TopologyAttr("resid", [1])

    # Pre-allocate optional arrays on the Timestep so assignment works.
    if _has_vel:
        u.trajectory.ts.has_velocities = True
        u.trajectory.ts._velocities = np.zeros((n_atoms, 3), dtype=np.float32)
    if _has_frc:
        u.trajectory.ts.has_forces = True
        u.trajectory.ts._forces = np.zeros((n_atoms, 3), dtype=np.float32)

    print(f"Writing TRR: {output_trr}")
    with TRRWriter(output_trr, n_atoms=n_atoms) as writer:
        for frame_idx, atoms in enumerate(frames):
            sim_time = frame_idx * skip * dt  # ps

            # ---- positions (Å  →  Å, MDAnalysis works in Å internally) ----
            positions = atoms.get_positions()  # shape (N, 3), Å

            # ---- velocities (ASE: Å/fs → MDAnalysis: Å/ps) ----------------
            velocities = None
            if _has_vel and atoms.get_velocities() is not None:
                velocities = atoms.get_velocities() * 1000.0  # Å/fs → Å/ps

            # ---- forces (ASE: eV/Å → MDAnalysis: kJ/(mol·Å)) --------------
            forces = None
            if _has_frc:
                try:
                    raw_forces = atoms.get_forces(apply_constraint=False)
                    # 1 eV/Å = 96.4853 kJ/(mol·Å)
                    forces = raw_forces * 96.4853
                except Exception:
                    pass

            # ---- periodic box (ASE cell in Å, orthorhombic or triclinic) --
            cell = atoms.get_cell()
            if cell.orthorhombic:
                lengths = cell.lengths()       # [a, b, c] in Å
                angles  = np.array([90.0, 90.0, 90.0])
            else:
                lengths = cell.lengths()
                angles  = cell.angles()        # [α, β, γ] in degrees

            # ---- update Universe positions for this frame ------------------
            u.atoms.positions = positions
            if velocities is not None:
                u.trajectory.ts._velocities[:] = velocities
            if forces is not None:
                u.trajectory.ts._forces[:] = forces

            # ---- set box dimensions [a, b, c, α, β, γ] --------------------
            u.trajectory.ts.dimensions = np.concatenate([lengths, angles])

            # ---- assign time & step ----------------------------------------
            u.trajectory.ts.time  = sim_time
            u.trajectory.ts.frame = frame_idx

            writer.write(u.atoms)

            if (frame_idx + 1) % max(1, n_frames // 10) == 0 or frame_idx == n_frames - 1:
                print(f"  Written frame {frame_idx + 1:>6d}/{n_frames}  (t = {sim_time:.1f} ps)")

    print(f"\nDone. Output written to: {output_trr}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert an ASE .traj file to a GROMACS .trr file."
    )
    parser.add_argument("input",  help="Input ASE trajectory (.traj)")
    parser.add_argument("output", help="Output GROMACS trajectory (.trr)")
    parser.add_argument(
        "--dt",
        type=float,
        default=1.0,
        help="Timestep between saved frames in picoseconds (default: 1.0 ps)",
    )
    parser.add_argument(
        "--skip",
        type=int,
        default=1,
        help="Write every Nth frame (default: 1, i.e. all frames)",
    )
    args = parser.parse_args()

    traj_to_trr(
        input_traj=args.input,
        output_trr=args.output,
        dt=args.dt,
        skip=args.skip,
    )


if __name__ == "__main__":
    main()