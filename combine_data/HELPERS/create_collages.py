import os
import math
import argparse
from PIL import Image, ImageDraw, ImageFont
from multiprocessing import Pool

def create_collage(image_paths, output_path, img_per_row=None, show_labels=True, padding=0):
    """
    Creates a collage from a list of image paths and saves it.
    Tries to make the grid as square as possible if img_per_row is not specified.
    Adds filenames as labels below each image if show_labels is True.
    """
    if not image_paths:
        print(f"No images found for {output_path}")
        return

    images = []
    for p in image_paths:
        try:
            img = Image.open(p)
            images.append(img)
        except Exception as e:
            print(f"Could not open image {p}: {e}")
            continue

    if not images:
        print(f"No valid images to process for {output_path}")
        return

    # Determine grid size
    num_images = len(images)
    if img_per_row is None:
        cols = math.ceil(math.sqrt(num_images))
    else:
        cols = img_per_row
    rows = math.ceil(num_images / cols)

    # Get max width and height for sizing, and add space for labels
    label_height = 30 if show_labels else 0  # Height for the filename label
    widths, heights = zip(*(i.size for i in images))
    max_width = max(widths)
    max_height = max(heights)
    
    cell_width = max_width + padding
    cell_height = max_height + label_height + padding

    collage_width = cols * cell_width
    collage_height = rows * cell_height

    collage_image = Image.new('RGB', (collage_width, collage_height), color='white')
    draw = ImageDraw.Draw(collage_image)
    
    try:
        # Use a common font if available, otherwise default
        font = ImageFont.truetype("arial.ttf", 15)
    except IOError:
        font = ImageFont.load_default()

    current_x, current_y = 0, 0
    for i, img in enumerate(images):
        # Calculate position to paste the image (centered in its cell if not max size)
        paste_x = current_x + (max_width - img.width) // 2
        paste_y = current_y + (max_height - img.height) // 2
        
        collage_image.paste(img, (paste_x, paste_y))

        if show_labels:
            # Add filename label
            filename = os.path.basename(image_paths[i])
            text_width = draw.textlength(filename, font=font)
            text_x = current_x + (cell_width - text_width) / 2
            text_y = current_y + max_height + 5 # 5 pixels padding
            draw.text((text_x, text_y), filename, fill="black", font=font)
        
        current_x += cell_width
        if (i + 1) % cols == 0:
            current_x = 0
            current_y += cell_height
            
    try:
        collage_image.save(output_path)
        print(f"Collage saved to {output_path}")
    except Exception as e:
        print(f"Error saving collage {output_path}: {e}")

def find_images_and_create_collages_precise(base_dir, show_labels=True, padding=0, exclude_list=None):
    """
    Finds DEGs*/plots/subdir and creates collages of .png images in each subdir. (Precise mode)
    """
    if exclude_list is None:
        exclude_list = []
    degs_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and d.startswith("DEGs")]

    for deg_dir_name in degs_dirs:
        deg_path = os.path.join(base_dir, deg_dir_name)
        plots_path = os.path.join(deg_path, "plots")

        if os.path.isdir(plots_path):
            print(f"Processing plots in: {plots_path}")
            for sub_plot_dir_name in os.listdir(plots_path):
                sub_plot_dir_path = os.path.join(plots_path, sub_plot_dir_name)
                if os.path.isdir(sub_plot_dir_path):
                    print(f"  Processing sub-directory: {sub_plot_dir_path}")
                    png_files = [os.path.join(sub_plot_dir_path, f) 
                                 for f in os.listdir(sub_plot_dir_path) 
                                 if f.lower().endswith(".png") and f.lower() != "collage.png" and f not in exclude_list]
                    
                    if png_files:
                        collage_output_path = os.path.join(sub_plot_dir_path, "collage.png")
                        create_collage(png_files, collage_output_path, show_labels=show_labels, padding=padding) 
                    else:
                        print(f"    No .png files found in {sub_plot_dir_path} (excluding 'collage.png')")
        else:
            print(f"No 'plots' directory found in {deg_path}")

def find_images_and_create_collages_general(data_dir, show_labels=True, padding=0, exclude_list=None):
    """
    Finds all .png files in data_dir and its subdirectories, creates a single collage,
    and saves it in data_dir. (General mode)
    """
    if exclude_list is None:
        exclude_list = []
    print(f"Processing all .png files in: {data_dir} (General mode)")
    png_files = []
    for root, _, files in os.walk(data_dir):
        for file in files:
            if file.lower().endswith(".png") and file.lower() != "collage.png" and file not in exclude_list:
                png_files.append(os.path.join(root, file))
    
    if png_files:
        collage_output_path = os.path.join(data_dir, "collage.png")
        create_collage(png_files, collage_output_path, show_labels=show_labels, padding=padding)
    else:
        print(f"No .png files found in {data_dir} (excluding 'collage.png')")

def process_directory(args_tuple):
    """Helper function to unpack arguments for multiprocessing."""
    data_dir, mode, show_labels, padding, exclude_list = args_tuple
    print(f"Running in '{mode}' mode on directory: {data_dir}")
    if mode == 'precise':
        find_images_and_create_collages_precise(data_dir, show_labels=show_labels, padding=padding, exclude_list=exclude_list)
    elif mode == 'general':
        find_images_and_create_collages_general(data_dir, show_labels=show_labels, padding=padding, exclude_list=exclude_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create collages from images in one or more directories in parallel. Can operate in two modes:\n"
                    "'precise': Looks for images in specific subdirectories (DEGs*/plots/...). \n"
                    "'general': Scans all subdirectories for .png files to create a single collage.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("data_dirs", nargs='+', help="One or more base directories to search for images.")
    parser.add_argument("--mode", choices=["precise", "general"], default="precise", 
                        help="Mode of operation:\n"
                             "  precise: Looks for DEGs*/plots/subdir patterns.\n"
                             "  general: Uses all .png files in data_dir recursively.")
    parser.add_argument("--no-labels", dest='show_labels', action="store_false",
                        help="Disable printing filenames as labels below images.")
    parser.add_argument("--padding", type=int, default=0,
                        help="Set padding in pixels between images (default: 0).")
    parser.add_argument("--exclude", nargs='+', default=[],
                        help="A list of filenames to exclude from the collage.")
    
    args = parser.parse_args()

    print("Starting collage creation process...")
    print("Make sure you have Pillow installed: pip install Pillow")

    # Prepare arguments for parallel processing
    process_args = [(data_dir, args.mode, args.show_labels, args.padding, args.exclude) for data_dir in args.data_dirs]

    # Use a Pool of workers to process directories in parallel
    with Pool() as pool:
        pool.map(process_directory, process_args)
    

    print("Collage creation process finished.")

    """
    You can now run the script from your terminal and provide a space-separated list of directories.

    ### Example Usage:

    **To run in `general` mode for multiple directories:**
    ```bash
    python create_collages.py --mode general "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_conditions/both_genotypes/Combined_GC_Mutant_vs_Control" "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_conditions/both_genotypes/Mature_GC_Mutant_vs_Control"
    ```

    **To run in `precise` mode for multiple directories:**
    ```bash
    python create_collages.py --mode precise /path/to/data1 /path/to/data2 /path/to/data3
    ```
    """

"""
python create_collages.py --mode general --no-labels \
    --exclude "GSEA_Focus_Neuro_CellDeath.png" -- \
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_conditions/both_genotypes/Combined_GC_Mutant_vs_Control" \
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_conditions/both_genotypes/Mature_GC_Mutant_vs_Control" \
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_conditions/emx1/Combined_GC_Mutant_vs_Control" \
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_conditions/emx1/Mature_GC_Mutant_vs_Control" \
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_conditions/nestin/Combined_GC_Mutant_vs_Control" \
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_conditions/nestin/Mature_GC_Mutant_vs_Control"
"""


"""
python create_collages.py --mode general --no-labels \
    --exclude "GSEA_Focus_Neuro_CellDeath.png" -- \
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/ctrl/both_genotypes/Combined_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/ctrl/both_genotypes/Mature_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/ctrl/emx1/Combined_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/ctrl/emx1/Mature_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/ctrl/nestin/Combined_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/ctrl/nestin/Mature_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/mut/both_genotypes/Combined_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/mut/both_genotypes/Mature_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/mut/emx1/Combined_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/mut/emx1/Mature_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/mut/nestin/Combined_GC"
    "D:/Github/SRF_Linda_RNA/combine_data/results_from_raw/gsea_analysis_between_clusters/mut/nestin/Mature_GC"
"""


