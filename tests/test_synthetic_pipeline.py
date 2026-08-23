from pathlib import Path
import subprocess
import sys


def test_synthetic_pipeline(tmp_path: Path):
    root = Path(__file__).resolve().parents[1]
    subprocess.run([
        sys.executable,
        str(root / "paper_repro" / "synthetic" / "generate_synthetic_example.py"),
        "--config", str(root / "paper_repro" / "synthetic" / "example_config.yaml"),
        "--output", str(tmp_path / "example"),
    ], cwd=root, check=True)
    assert (tmp_path / "example" / "synthetic_bounds.csv").is_file()
    assert (tmp_path / "example" / "synthetic_bounds.pdf").is_file()
