# SpatialFilter


This code is for filtering cells in spatial transcriptomics datasets using binary mask images. 


## Example

```
# The mask must have white regions for the cells you want and black background. The resolution should be exactly the same as the original tissue_lowres_image.png
result <- filter_cells_by_mask(spatial.obj, "./spatial/tissue_lowres_image_mask.png")

# Get the filtered Seurat object
filtered_obj <- result$filtered_seurat
```


Overall this function returns a list with these key components:

filtered_seurat: The Seurat object with only cells inside the mask
in_mask: Logical vector indicating which cells are inside the mask
cells_in_mask: Cell IDs of all cells inside the mask
cells_outside_mask: Cell IDs of all cells outside the mask
mask: The binary mask image used for filtering
scaling: Alignment parameters (x_scale, y_scale, x_offset, y_offset)
stats: Summary statistics (total_cells, cells_in_mask, cells_outside_mask, percent_in_mask)
plot: Visualization of cells by mask status (if visualize=TRUE)
