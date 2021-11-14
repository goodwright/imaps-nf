import pandas as pd
import argparse


def annotation_to_barcode(annotation_path, output_path="barcode.csv"):
    """
    Extracts the barcode from the annotation
    input: annotation file as a string path
    """
    data = pd.read_csv(annotation_path)
    data = data[['Sample Name', "5' Barcode"]]
    data.to_csv(output_path, index=False, header=False)


def main():
    parser = argparse.ArgumentParser(description="Extracts the barcode from the annotation")
    parser.add_argument("-a", "--annotation", help="Path to the annotation file", required=True)
    parser.add_argument("-o", "--output", help="Path to the output file", required=False)
    args = parser.parse_args()
    annotation_to_barcode(args.annotation)

if __name__ == "__main__":
    main()