import argparse
def parse_args():
    parser = argparse.ArgumentParser(description = "Added Mass")
    parser.add_argument("--a", type = float, default = 1, help = "Major axis of ellipse, half length of square")
    parser.add_argument("--b", type = float, default = 1, help = "Minor axis of ellipse")
    parser.add_argument("--N", type = int, default = 32, help = "Number of boundary nodes")
    parser.add_argument("--number", type = int, default = 5, help = "Number of times to evaluate added mass")
    return parser.parse_args()
