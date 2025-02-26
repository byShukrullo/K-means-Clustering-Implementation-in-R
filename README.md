#' @param data A numeric matrix or data frame where rows are observations and columns are variables
#' @param k The number of clusters to form
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @return A list containing cluster assignments and centroids
kmeans_custom <- function(data, k, max_iter = 100, tol = 1e-4) {
  # Convert data to matrix if it's a data frame
  if(is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  # Get dimensions
  n <- nrow(data)
  p <- ncol(data)
  
  # Randomly initialize k centroids
  set.seed(42)  # For reproducibility
  centroid_indices <- sample(1:n, k)
  centroids <- data[centroid_indices, , drop = FALSE]
  
  # Initialize cluster assignments
  cluster <- rep(0, n)
  
  # Main K-means loop
  for(iter in 1:max_iter) {
    # Store previous centroids to check for convergence
    old_centroids <- centroids
    
    # Assign each point to nearest centroid
    for(i in 1:n) {
      # Calculate distance to each centroid
      distances <- apply(centroids, 1, function(cent) sum((data[i, ] - cent)^2))
      # Assign to closest centroid
      cluster[i] <- which.min(distances)
    }
    
    # Update centroids
    for(j in 1:k) {
      if(sum(cluster == j) > 0) {
        centroids[j, ] <- colMeans(data[cluster == j, , drop = FALSE])
      }
    }
    
    # Check for convergence
    if(max(colSums((centroids - old_centroids)^2)) < tol) {
      break
    }
  }
  
  # Calculate within-cluster sum of squares
  within_ss <- numeric(k)
  for(j in 1:k) {
    if(sum(cluster == j) > 0) {
      within_ss[j] <- sum(apply(data[cluster == j, , drop = FALSE], 1, function(x) {
        sum((x - centroids[j, ])^2)
      }))
    }
  }
  
  # Calculate total within-cluster sum of squares
  total_within_ss <- sum(within_ss)
  
  # Return results
  return(list(
    cluster = cluster,
    centroids = centroids,
    within_cluster_ss = within_ss,
    total_within_ss = total_within_ss,
    iterations = iter
  ))
}

#' Find optimal k using elbow method
#' 
#' @param data A numeric matrix or data frame
#' @param max_k Maximum number of clusters to try
#' @return A plot showing total within-cluster SS for different k values
find_optimal_k <- function(data, max_k = 10) {
  wss <- numeric(max_k)
  
  for(k in 1:max_k) {
    km <- kmeans_custom(data, k)
    wss[k] <- km$total_within_ss
  }
  
  # Plot elbow curve
  plot(1:max_k, wss, type = "b", pch = 19,
       xlab = "Number of clusters (k)",
       ylab = "Total within-cluster sum of squares",
       main = "Elbow Method for Optimal k")
  
  return(wss)
}

#' Plot K-means clusters in 2D
#' 
#' @param data A numeric matrix or data frame with 2 columns
#' @param clusters Cluster assignments
#' @param centroids Matrix of cluster centroids
plot_kmeans_clusters <- function(data, clusters, centroids) {
  if(ncol(data) != 2) {
    stop("Data must have exactly 2 columns for 2D plotting")
  }
  
  # Create a color palette for the clusters
  colors <- rainbow(nrow(centroids))
  
  # Plot the data points
  plot(data, col = colors[clusters], pch = 20, cex = 1.5,
       main = "K-means Clustering Result",
       xlab = "Feature 1", ylab = "Feature 2")
  
  # Add centroids
  points(centroids, col = "black", pch = 8, cex = 2, lwd = 2)
  
  # Add a legend
  legend("topright", legend = paste("Cluster", 1:nrow(centroids)), 
         col = colors, pch = 20, cex = 0.8)
}

# Example 1: Using the custom K-means implementation
# Generate example data
set.seed(123)
example_data <- rbind(
  matrix(rnorm(100, mean = 0, sd = 0.3), ncol = 2),
  matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2),
  matrix(rnorm(100, mean = 2, sd = 0.3), ncol = 2)
)
colnames(example_data) <- c("x", "y")

# Run custom K-means
k <- 3
km_result <- kmeans_custom(example_data, k)

# Plot the results
plot_kmeans_clusters(example_data, km_result$cluster, km_result$centroids)

# Example 2: Using built-in kmeans function in R
# For comparison with the custom implementation
kmeans_builtin <- kmeans(example_data, centers = k, nstart = 25)

# Print results comparison
cat("Custom K-means total within-cluster SS:", km_result$total_within_ss, "\n")
cat("Built-in K-means total within-cluster SS:", kmeans_builtin$tot.withinss, "\n")

# Example 3: Find optimal k using elbow method
wss <- find_optimal_k(example_data, max_k = 10)

# Example 4: Working with real-world dataset - Iris
data(iris)
iris_data <- iris[, 1:4]  # Extract numeric columns

# Find optimal k for Iris
wss_iris <- find_optimal_k(iris_data, max_k = 10)

# Apply K-means with k=3 (known number of species)
km_iris <- kmeans_custom(iris_data, k = 3)

# Evaluate clustering against actual species
table(km_iris$cluster, iris$Species)

# Example 5: PCA visualization for higher dimensional data
iris_pca <- prcomp(iris_data, scale. = TRUE)
iris_pca_data <- as.data.frame(iris_pca$x[, 1:2])  # First two principal components

# Run K-means on PCA data
km_iris_pca <- kmeans_custom(iris_pca_data, k = 3)

# Plot results with actual species information
plot_kmeans_pca <- function() {
  # Set up the plot
  plot(iris_pca_data$PC1, iris_pca_data$PC2, 
       col = km_iris_pca$cluster, pch = as.numeric(iris$Species),
       main = "Iris PCA with K-means Clusters",
       xlab = "PC1", ylab = "PC2")
  
  # Add centroids
  points(km_iris_pca$centroids, col = "black", pch = 8, cex = 2, lwd = 2)
  
  # Add a legend for clusters and species
  legend("topright", 
         legend = c(paste("Cluster", 1:3), levels(iris$Species)),
         col = c(1:3, 1), 
         pch = c(1, 1, 1, 1:3),
         cex = 0.8)
}

plot_kmeans_pca()

# Example 6: Silhouette analysis for cluster quality
library(cluster)

silhouette_kmeans <- function(data, k) {
  # Run K-means
  km <- kmeans(data, centers = k, nstart = 25)
  
  # Calculate silhouette
  sil <- silhouette(km$cluster, dist(data))
  
  # Plot silhouette
  plot(sil, main = paste("Silhouette Plot for k =", k))
  
  # Return average silhouette width
  return(mean(sil[, 3]))
}

# Compare silhouette scores for different k values
silhouette_scores <- sapply(2:6, function(k) silhouette_kmeans(iris_data, k))
plot(2:6, silhouette_scores, type = "b", xlab = "Number of clusters (k)", 
     ylab = "Average Silhouette Width", main = "Silhouette Analysis")

# Example 7: K-means with multiple random starts
kmeans_multiple_starts <- function(data, k, n_starts = 10) {
  best_wss <- Inf
  best_result <- NULL
  
  for(i in 1:n_starts) {
    result <- kmeans_custom(data, k)
    if(result$total_within_ss < best_wss) {
      best_wss <- result$total_within_ss
      best_result <- result
    }
  }
  
  return(best_result)
}

km_multiple <- kmeans_multiple_starts(example_data, k = 3, n_starts = 10)
plot_kmeans_clusters(example_data, km_multiple$cluster, km_multiple$centroids)
