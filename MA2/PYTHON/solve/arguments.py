import argparse
def parse_args():
    parser = argparse.ArgumentParser(description = "Added Mass")
    parser.add_argument("--Nx", type = int, default = 100, help = "Number of points in x")
    parser.add_argument("--Ny", type = int, default = 25, help = "Number of points in y")
    parser.add_argument("--L", type = float, default = 2, help = "Length, along x")
    parser.add_argument("--D", type = float, default = 1, help = "Draught, along y")
    parser.add_argument("--kD", type = float, default = 1.2, help = "Wavenumber times draught")
    return parser.parse_args()