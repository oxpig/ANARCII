# configuration_utils.py
import os

import psutil
import torch


def configure_cpus(ncpu: int) -> int:
    """Configure CPU allocation based on availability and user input."""
    available_cpus = get_available_cpus()
    if ncpu == -1:
        torch.set_num_threads(available_cpus)
        return available_cpus
    elif ncpu <= available_cpus:
        torch.set_num_threads(ncpu)
        return ncpu
    else:
        print("Number of requested CPUs exceeds available CPUs.")
        torch.set_num_threads(available_cpus)
        return available_cpus


def configure_device(cpu: bool, ncpu: int) -> torch.device:
    """
    Configure computation device (CPU or GPU).
    It will default to CPU if no GPU is available or if cpu=True.
    """
    if cpu or not torch.cuda.is_available():
        device = torch.device("cpu")
    else:
        device = torch.device("cuda")

    print(f"Using device {str(device).upper()} with {ncpu} CPUs")
    return device


def get_available_cpus():
    try:
        process = psutil.Process()
        cpu_affinity = process.cpu_affinity()
        return len(cpu_affinity)
    except Exception:  # noqa: BLE001
        if "SLURM_CPUS_PER_TASK" in os.environ:
            return int(os.environ["SLURM_CPUS_PER_TASK"])
        elif "PBS_NP" in os.environ:
            return int(os.environ["PBS_NP"])
        elif "LSB_DJOB_NUMPROC" in os.environ:
            return int(os.environ["LSB_DJOB_NUMPROC"])
        # Final fallback: detect all cores
        return os.cpu_count()
