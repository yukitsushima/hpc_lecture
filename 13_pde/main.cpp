#include <cstdio>
#include <vector>
#include <tuple>
#include <string>
#include <fstream>

//For debug
void print_2dvec(std::vector<std::vector<double>> a,int nx, int ny) {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny-1; j++) {
      printf("%f,", a[i][j]);
    }
    printf("%f\n", a[i][ny-1]);
  }
}

//For release
void write_2dvec(std::vector<std::vector<double>> a,int nx, int ny,std::string filename) {
  std::ofstream writing_file(filename);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny-1; j++) {
      writing_file<<a[i][j]<<",";
    }
    writing_file<<a[i][ny-1]<<"\n";
  }
  writing_file.close();
}

std::vector<std::vector<double>> build_up_b(std::vector<std::vector<double>> b,double rho,double dt,std::vector<std::vector<double>> u,std::vector<std::vector<double>> v,double dx,double dy,int nx,int ny){
  std::vector<std::vector<double>> b_next(ny,std::vector<double>(nx,0));
  for (int x = 1; x < nx-1; x++) {
    for (int y = 1; y < ny-1; y++) {
      b_next[x][y] =
      (
      (((u[x+1][y]-u[x-1][y])/(2*dx))+((v[x][y+1]-v[x][y-1])/(2*dy)))/dt
      -(((u[x+1][y]-u[x-1][y])/(2*dx))*((u[x+1][y]-u[x-1][y])/(2*dx)))
      -2*(((u[x][y+1]-u[x][y-1])/(2*dy))*((v[x+1][y]-v[x-1][y])/(2*dx)))
      -(((v[x][y+1]-v[x][y-1])/(2*dy))*((v[x][y+1]-v[x][y-1])/(2*dy)))
      )*rho;
    }
  }
  return b_next;
}
std::vector<std::vector<double>> pressure_poisson(std::vector<std::vector<double>> p,double dx,double dy,std::vector<std::vector<double>> b,int nx,int ny,int nit){
  std::vector<std::vector<double>> p_next(ny,std::vector<double>(nx,0));
  for (int ctr = 0; ctr < nit; ctr++) {
    for (int x = 1; x < nx-1; x++) {
      for (int y = 1; y < ny-1; y++) {
        p_next[x][y] =
        ((p[x+1][y]+p[x-1][y])*(dy*dy)+(p[x][y+1]+p[x][y-1])*(dx*dx))/(2*(dx*dx+dy*dy))
        -(dx*dx*dy*dy)/(2*(dx*dx+dy*dy))*b[x][y];
      }
    }
    for (int x = 0; x < nx; x++) {
      p_next[x][ny-1] = p_next[x][ny-2];
      p_next[x][0] = p_next[x][1];
    }
    for (int y = 0; y < ny; y++) {
      p_next[0][y] = p_next[1][y];
      p_next[nx-1][y] = 0;
    }
    p = p_next;
  }
  return p_next;
}

std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>>
cavity_flow(int nt,std::vector<std::vector<double>> u,std::vector<std::vector<double>> v,
  double dt,double dx,double dy, std::vector<std::vector<double>> p,double rho,double nu,int nx, int ny,int nit) {
  std::vector<std::vector<double>> b(ny,std::vector<double>(nx,0));
  for (int i = 0; i < nt; i++) {
    std::vector<std::vector<double>> u_next(ny,std::vector<double>(nx,0));
    std::vector<std::vector<double>> v_next(ny,std::vector<double>(nx,0));
    b = build_up_b(b,rho,dt,u,v,dx,dy,nx,ny);
    p = pressure_poisson(p,dx,dy,b,nx,ny,nit);
    //Process u,v,p -> u_next,v_next,p_next
    for (int x = 1; x < nx-1; x++) {
      for (int y = 1; y < ny-1; y++) {
        u_next[x][y] =
        u[x][y]-u[x][y]*(dt/dx)*(u[x][y]-u[x-1][y])-v[x][y]*(dt/dy)*(u[x][y]-u[x][y-1])
        -(dt/(rho*2*dx))*(p[x+1][y]-p[x-1][y])
        +nu*((dt/(dx*dx))*(u[x+1][y]-2*u[x][y]+u[x-1][y])+(dt/(dy*dy))*(u[x][y+1]-2*u[x][y]+u[x][y-1]));
        v_next[x][y] = v[x][y]-u[x][y]*(dt/dx)*(v[x][y]-v[x-1][y])-v[x][y]*(dt/dy)*(v[x][y]-v[x][y-1])
        -(dt/(2*rho*dy))*(p[x][y+1]-p[x][y-1])
        +nu*((dt/(dx*dx))*(v[x+1][y]-2*v[x][y]+v[x-1][y])+(dt/(dy*dy))*(v[x][y+1]-2*v[x][y]+v[x][y-1]));
      }
    }
    for (int y = 0; y < ny; y++) {
      u_next[0][y] = 0;
      v_next[0][y] = 0;
      u_next[nx-1][y] = 1;
      v_next[nx-1][y] = 0;
    }
    for (int x = 0; x < nx; x++) {
      u_next[x][0] = 0;
      u_next[x][ny-1] = 0;
      v_next[x][0] = 0;
      v_next[x][ny-1] = 0;
    }
    u = u_next;
    v = v_next;
  }
  return std::forward_as_tuple(u,v,p);
}

int main(int argc, char** argv) {
  const int nx = 41;
  const int ny = 41;
  const int nt = 100;
  const int nit = 50;
  const int c = 1;
  const double dx = 2.0/(nx-1);
  const double dy = 2.0/(ny-1);
  const double rho = 1;
  const double nu = 0.1;
  const double dt = 0.001;
  std::vector<std::vector<double>> u(ny,std::vector<double>(nx,0));
  std::vector<std::vector<double>> v(ny,std::vector<double>(nx,0));
  std::vector<std::vector<double>> p(ny,std::vector<double>(nx,0));
  std::tie(u,v,p) = cavity_flow(nt,u,v,dt,dx,dy,p,rho,nu,nx,ny,nit);
  //print_2dvec(u,nx,ny);
  //print_2dvec(v,nx,ny);
  write_2dvec(u,nx,ny,"u.csv");
  write_2dvec(v,nx,ny,"v.csv");
}
