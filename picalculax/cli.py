"""Command-line interface for pICalculax."""

import argparse
from .core import pred_fasta, pred_pln

def main():
    """Main entry point for the command-line interface."""
    parser = argparse.ArgumentParser(
        description='Predict isoelectric point pI of peptides and modified peptides'
    )
    parser.add_argument(
        '--fasta',
        nargs='+',
        type=str,
        help='Predict fasta sequence'
    )
    parser.add_argument(
        '--pln',
        nargs='+',
        type=str,
        help='Predict PLN sequence (Requires Proteax Desktop)'
    )
    args = parser.parse_args()

    # Predict lists
    if args.fasta is not None:
        pred_fasta(args.fasta)
    if args.pln is not None:
        pred_pln(args.pln)

if __name__ == '__main__':
    main() 