import csv
import re

def clean_csv_column(input_filepath, output_filepath):
    """
    Reads a CSV file, removes leading numbers and a space from the second column,
    and writes the cleaned data to a new CSV file.
    """
    with open(input_filepath, 'r', newline='') as infile, \
         open(output_filepath, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        # Write the header row as is
        header = next(reader)
        writer.writerow(header)

        # Process the rest of the rows
        for row in reader:
            if len(row) > 1:
                # Remove leading numbers and the following space from the second column
                row[1] = re.sub(r'^\d+\s*', '', row[1])
            writer.writerow(row)

if __name__ == "__main__":
    input_file = '../combine_data/models/second_layer_LB.csv'
    output_file = '../combine_data/models/second_layer_LB_cleaned.csv'
    clean_csv_column(input_file, output_file)
    print(f"Cleaned data written to {output_file}")

    input_file = '../combine_data/models/first_layer_LB.csv'
    output_file = '../combine_data/models/first_layer_LB_cleaned.csv'
    clean_csv_column(input_file, output_file)
    print(f"Cleaned data written to {output_file}")