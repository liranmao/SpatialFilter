# SpatialFilter


This code is for filtering cells in spatial transcriptomics datasets using binary mask images. 


## How to use

```
# The mask must have white regions for the cells you want and black background. The resolution should be exactly the same as the original tissue_lowres_image.png
result <- filter_cells_by_mask(spatial.obj, "./spatial/tissue_lowres_image_mask.png")

# Get the filtered Seurat object
filtered_obj <- result$filtered_seurat
```

## Outputs
Overall this function returns a list with these key components:

- filtered_seurat: The Seurat object with only cells inside the mask
- in_mask: Logical vector indicating which cells are inside the mask
- cells_in_mask: Cell IDs of all cells inside the mask
- cells_outside_mask: Cell IDs of all cells outside the mask
- mask: The binary mask image used for filtering
- scaling: Alignment parameters (x_scale, y_scale, x_offset, y_offset)
- stats: Summary statistics (total_cells, cells_in_mask, cells_outside_mask, percent_in_mask)
- plot: Visualization of cells by mask status (if visualize=TRUE)


## Example

We are interested in the spine development and only want cells from that region. Ther following is the image from spatial folder in standard analysis pipeline:

<img width="200" alt="image" src="https://github.com/user-attachments/assets/8dc91ab8-aae4-4063-a9b0-a7ed101d9897" />

Use PS to draw mask manully:

<img width="200" alt="image" src="https://github.com/user-attachments/assets/70e9ff08-d1fe-4e55-9e80-f308c7374ac5" />

Before filtering,

<img width="220" alt="image" src="https://github.com/user-attachments/assets/8d9a4999-e204-4780-b424-bcb04f6a9181" />


Using filter_cells_by_mask function, you can get cells:

<img width="220" alt="image" src="https://github.com/user-attachments/assets/6cc74841-010a-4273-97f0-713bfa9faa5c" />


