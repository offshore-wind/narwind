#ifndef GEODESIC_BFS
#define GEODESIC_BFS
#include <queue>
#include <iostream>
#include <string>

// [[Rcpp::export]]
std::int64_t create_hash (double x0, double y0, double x1, double y1){
  auto d0 = (x0 * y0) + std::floor(std::pow((std::abs(x0 - y0) - 1), 2))/4;
  auto d1 = (x1 * y1) + std::floor(std::pow((std::abs(x1 - y1) - 1), 2))/4;
  return std::trunc(d0 * d1);
}

// std::vector<double> cellFromXY (std::vector<double> x, std::vector<double> y, double missing) {
//   // size of x and y should be the same
//   
//   size_t size = x.size();
//   std::vector<double> cells(size);
//   
//   SpatExtent extent = getExtent();
//   double yr_inv = nrow() / (extent.ymax - extent.ymin);
//   double xr_inv = ncol() / (extent.xmax - extent.xmin);
//   
//   for (size_t i = 0; i < size; i++) {
//     // cannot use trunc here because trunc(-0.1) == 0
//     long row = std::floor((extent.ymax - y[i]) * yr_inv);
//     // points in between rows go to the row below
//     // except for the last row, when they must go up
//     if (y[i] == extent.ymin) {
//       row = nrow()-1 ;
//     }
//     
//     long col = std::floor((x[i] - extent.xmin) * xr_inv);
//     // as for rows above. Go right, except for last column
//     if (x[i] == extent.xmax) {
//       col = ncol() - 1 ;
//     }
//     long nr = nrow();
//     long nc = ncol();
//     if (row < 0 || row >= nr || col < 0 || col >= nc) {
//       cells[i] = missing;
//     } else {
//       cells[i] = row * ncol() + col;
//     }
//   }
//   
//   return cells;
// }


// QItem for current location and distance 
// from source location 
class QItem { 
public: 
  int row; 
  int col; 
  int dist; 
  QItem(int x, int y, int w) 
    : row(x), col(y), dist(w) 
  { 
  } 
};

// [[Rcpp::export]]
int geoD(Eigen::MatrixXd mat,
         double x0,
         double y0,
         double x1,
         double y1,
         Eigen::VectorXd limits,
         Eigen::VectorXd resolution)
{
  // Retrieve row and col of cell corresponding to x,y
  std::size_t c0 = std::round((x0 - limits[0])/resolution[0]);
  std::size_t r0 = std::round((limits[3] - y0)/resolution[1]);
  std::size_t c1 = std::round((x1 - limits[0])/resolution[0]);
  std::size_t r1 = std::round((limits[3] - y1)/resolution[1]);
  
  // std::cout << r0 << " " << c0 << "\n";
  // std::cout << r1 << " " << c1 << "\n";
  
  // Source on land
  if (mat(r0,c0) == 0){
    // std::cout << "Source on land!\n";
    return -1;
  }
  
  mat(r0,c0) = 2; // Source
  mat(r1,c1) = 3; // Destination
  
  // std::cout << mat(r1,c1) << "\n";
  
  int nr = mat.rows();
  int nc = mat.cols();
  
  QItem source(r0,c0,0); 
  
  // To keep track of visited QItems. 
  // Marking blocked cells as visited. 
  bool visited[nr][nc]; 
  for (int i = 0; i < nr; i++) { 
    for (int j = 0; j < nc; j++) 
    { 
      if (mat(i,j) == 0) 
        visited[i][j] = true; 
      else
        visited[i][j] = false;
    } 
    
  } 
  
  // Applying Breadth-First-Search (BFS) on matrix cells starting from source 
  std::queue<QItem> q;
  q.push(source); // Insert a new element at the end of the queue, after its current last element.
  visited[source.row][source.col] = true;
  
  // While size of queue > 0
  while (!q.empty()) {
    
    QItem p = q.front(); // Return reference to next element
    q.pop(); // Remove the next element in the queue, effectively reducing its size by one
    
    // Destination found;
    if (mat(p.row, p.col) == 3){
      // std::cout << "Destination found!\n";
      return p.dist * resolution[0];
    }

    // moving up
    if (p.row - 1 >= 0 && visited[p.row - 1][p.col] == false) {
      q.push(QItem(p.row - 1, p.col, p.dist + 1));
      visited[p.row - 1][p.col] = true;
    }

    
    // moving down
    if (p.row + 1 < nr && visited[p.row + 1][p.col] == false) {
      q.push(QItem(p.row + 1, p.col, p.dist + 1));
      visited[p.row + 1][p.col] = true;
    }
    
    // moving left
    if (p.col - 1 >= 0 && visited[p.row][p.col - 1] == false) {
      q.push(QItem(p.row, p.col - 1, p.dist + 1));
      visited[p.row][p.col - 1] = true;
    }
    
    // moving right
    if (p.col + 1 < nc && visited[p.row][p.col + 1] == false) {
      q.push(QItem(p.row, p.col + 1, p.dist + 1));
      visited[p.row][p.col + 1] = true;
    }
  }
  
  return -1;
  
}

// [[Rcpp::export]]
int geoDist(Eigen::MatrixXd mat,
            double x0,
            double y0,
            double x1,
            double y1,
            Eigen::VectorXd limits,
            Eigen::VectorXd resolution,
            int r)
{
  // Retrieve row and col of cell corresponding to x,y
  std::size_t c0 = std::round((x0 - limits[0])/resolution[0]);
  std::size_t r0 = std::round((limits[3] - y0)/resolution[1]);
  std::size_t c1 = std::round((x1 - limits[0])/resolution[0]);
  std::size_t r1 = std::round((limits[3] - y1)/resolution[1]);
  
  // std::cout << r0 << " " << c0 << "\n";
  // std::cout << r1 << " " << c1 << "\n";
  
  // Source on land
  if (mat(r0,c0) == 0){
    // std::cout << "Source on land!\n";
    return -1;
  }
  
  mat(r0,c0) = 2; // Source
  mat(r1,c1) = 3; // Destination
  
  // std::cout << mat(r1,c1) << "\n";
  
  int nr = mat.rows();
  int nc = mat.cols();
  
  QItem source(r0,c0,0); 
  
  // To keep track of visited QItems. 
  // Marking blocked cells as visited. 
  bool visited[nr][nc];
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++)
    {
      if (mat(i,j) == 0)
        visited[i][j] = true;
      else
        visited[i][j] = false;
    }
  }

  // Applying Breadth-First-Search (BFS) on matrix cells starting from source 
  std::queue<QItem> q;
  q.push(source); // Insert a new element at the end of the queue, after its current last element.
  visited[source.row][source.col] = true;
  
  int i = 0;
  
  // While size of queue > 0
  while (!q.empty()) {
    
    // while(i < r) {
    if(i == r){
      // std::cout << "Stop!\n";
      break;
    }
    
    i++;
    QItem p = q.front(); // Return reference to next element
    q.pop(); // Remove the next element in the queue, effectively reducing its size by one
    
    // Destination found;
    if (mat(p.row, p.col) == 3){
      // std::cout << "Destination found!\n";
      return p.dist * resolution[0];
    }
    
    // moving up
    if (p.row - 1 >= 0 && visited[p.row - 1][p.col] == false) {
      q.push(QItem(p.row - 1, p.col, p.dist + 1));
      visited[p.row - 1][p.col] = true;
    }
    
    // moving down
    if (p.row + 1 < nr && visited[p.row + 1][p.col] == false) {
      q.push(QItem(p.row + 1, p.col, p.dist + 1));
      visited[p.row + 1][p.col] = true;
    }
    
    // moving left
    if (p.col - 1 >= 0 && visited[p.row][p.col - 1] == false) {
      q.push(QItem(p.row, p.col - 1, p.dist + 1));
      visited[p.row][p.col - 1] = true;
    }
    
    // moving right
    if (p.col + 1 < nc && visited[p.row][p.col + 1] == false) {
      q.push(QItem(p.row, p.col + 1, p.dist + 1));
      visited[p.row][p.col + 1] = true;
    }
    
  }
  
  return -1;
  
}

#endif