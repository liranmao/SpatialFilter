#' Filter cells in a Seurat spatial object based on a binary mask image
#'
#' @param seurat_obj A Seurat object with spatial data
#' @param mask_path Path to the binary mask image
#' @param use_image_coords Whether to use imagerow/imagecol (TRUE) or row/col (FALSE)
#' @param invert Whether to invert the mask (TRUE if white is background, FALSE if white is foreground)
#' @param x_offset Manual x-axis offset to align coordinates
#' @param y_offset Manual y-axis offset to align coordinates
#' @param x_scale Manual x-axis scaling factor (if NULL, auto-calculated)
#' @param y_scale Manual y-axis scaling factor (if NULL, auto-calculated)
#' @param visualize Whether to create diagnostic plots
#'
#' @return A list containing the filtered Seurat object, mask status, and visualization
#'
filter_cells_by_mask <- function(seurat_obj, 
                                 mask_path, 
                                 use_image_coords = TRUE,
                                 invert = FALSE,
                                 x_offset = 0,
                                 y_offset = 0,
                                 x_scale = NULL,
                                 y_scale = NULL,
                                 visualize = TRUE) {
  
  # Check required packages
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop("Package 'EBImage' is required. Install it with BiocManager::install('EBImage')")
  }
  
  # Load the mask image
  mask <- EBImage::readImage(mask_path)
  
  # Convert to grayscale if RGB
  if (length(dim(mask)) > 2) {
    mask <- EBImage::channel(mask, "gray")
  }
  
  # Ensure mask is binary (0=background, 1=foreground)
  if (max(mask) > 1) {
    mask <- mask > 0.5  # Convert to binary
  }
  
  # Invert if requested
  if (invert) {
    mask <- 1 - mask
  }
  
  # Get dimensions
  dim_mask <- dim(mask)
  
  # Get coordinates from Seurat object
  coords <- seurat_obj@images$slice1@coordinates
  
  if (use_image_coords) {
    # Use imagerow and imagecol
    coord_row <- coords$imagecol
    coord_col <- coords$imagerow
  } else {
    # Use row and col
    coord_row <- coords$col
    coord_col <- coords$row
  }
  
  # Calculate scaling factors if not provided
  if (is.null(x_scale)) {
    x_scale <- dim_mask[2] / max(coord_col)
  }
  
  if (is.null(y_scale)) {
    y_scale <- dim_mask[1] / max(coord_row)
  }
  
  # Scale coordinates to match mask dimensions and apply offsets
  mask_row <- round(coord_row * y_scale + y_offset)
  mask_col <- round(coord_col * x_scale + x_offset)
  
  # Ensure coordinates are within image boundaries
  valid_coords <- mask_row >= 1 & mask_row <= dim_mask[1] & 
    mask_col >= 1 & mask_col <= dim_mask[2]
  
  if (sum(!valid_coords) > 0) {
    warning(paste0(sum(!valid_coords), " cells have coordinates outside the mask boundaries and will be excluded"))
  }
  
  # Initialize mask status vector
  in_mask <- rep(FALSE, nrow(coords))
  
  # Check only valid coordinates
  valid_indices <- which(valid_coords)
  for (i in valid_indices) {
    in_mask[i] <- mask[mask_row[i], mask_col[i]] > 0
  }
  
  # Subset the Seurat object
  cells_in_mask <- rownames(coords)[in_mask]
  filtered_seurat <- subset(seurat_obj, cells = cells_in_mask)
  
  # Create result list
  result <- list(
    filtered_seurat = filtered_seurat,
    in_mask = in_mask,
    cells_in_mask = cells_in_mask,
    cells_outside_mask = rownames(coords)[!in_mask],
    mask = mask,
    scaling = list(x_scale = x_scale, y_scale = y_scale,
                   x_offset = x_offset, y_offset = y_offset)
  )
  
  # Add statistics
  result$stats <- list(
    total_cells = nrow(coords),
    cells_in_mask = sum(in_mask),
    cells_outside_mask = sum(!in_mask),
    percent_in_mask = round(100 * sum(in_mask) / nrow(coords), 2)
  )
  
  # Add visualization if requested
  if (visualize && requireNamespace("ggplot2", quietly = TRUE)) {
    # Create data frame for plotting
    plot_data <- data.frame(
      cell_id = rownames(coords),
      x = coord_row,
      y = coord_col,
      mask_x = mask_row,
      mask_y = mask_col,
      in_mask = factor(in_mask)
    )
    
    # Plot cells colored by mask status
    p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x=x, y=y, color=in_mask)) +
      ggplot2::geom_point(alpha=0.7) +
      ggplot2::scale_color_manual(values=c("TRUE"="blue", "FALSE"="red")) +
      ggplot2::labs(title = "Cells by mask status (original coordinates)", 
                    x = "Column", y = "Row") +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_minimal()
    
    result$plot <- p1
  }
  
  return(result)
}
