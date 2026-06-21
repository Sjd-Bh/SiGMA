import argparse
import os
import random

BASES = b"ATGC"

def write_random_fasta(
    filename,
    seq_name,
    length_bp,
    line_width=70,
    chunk_bp=1_000_000,  # 1 Mb chunks
):
    with open(filename, "wb") as f:
        f.write(f">{seq_name}\n".encode())

        remaining = length_bp
        buffer = b""

        while remaining > 0:
            n = min(chunk_bp, remaining)

            # Generate random bases
            bases = bytes(random.choices(BASES, k=n))
            buffer += bases
            remaining -= n

            # Flush complete lines
            while len(buffer) >= line_width:
                f.write(buffer[:line_width] + b"\n")
                buffer = buffer[line_width:]

        # Flush remainder
        if buffer:
            f.write(buffer + b"\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate random reference genome FASTA files efficiently."
    )
    parser.add_argument(
        "-length",
        type=int,
        nargs="+",
        required=True,
        help="Sequence lengths in kb (e.g. 200 400 600)"
    )
    parser.add_argument(
        "-o",
        required=True,
        help="Output directory"
    )
    args = parser.parse_args()

    os.makedirs(args.o, exist_ok=True)

    for length_kb in args.length:
        length_bp = length_kb * 1000
        out_fasta = os.path.join(
            args.o, f"reference_sequence_{length_kb}kb.fasta"
        )

        write_random_fasta(
            filename=out_fasta,
            seq_name=f"{length_kb}kb",
            length_bp=length_bp
        )

    print("FASTA generation completed successfully.")


if __name__ == "__main__":
    main()
