#include <cstdio>
#include <vector>
#include <tuple>
#include <string>
#include <fstream>
#include <iomanip>
#include <mpi.h>

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
  writing_file << std::fixed;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny-1; j++) {
      writing_file<<std::setprecision(15)<<a[i][j]<<",";
    }
    writing_file<<std::setprecision(15)<<a[i][ny-1]<<"\n";
  }
  writing_file.close();
}

//Same as Python code (transpose)
void write_2dvecT(std::vector<std::vector<double>> a,int nx, int ny,std::string filename) {
  std::ofstream writing_file(filename);
  writing_file << std::fixed;
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx-1; i++) {
      writing_file<<std::setprecision(15)<<a[i][j]<<",";
    }
    writing_file<<std::setprecision(15)<<a[nx-1][j]<<"\n";
  }
  writing_file.close();
}

std::vector<std::vector<double>>
build_up_b(double rho,double dt,std::vector<std::vector<double>> u,std::vector<std::vector<double>> v,double dx,double dy,int nx,int ny,int size, int rank){
  std::vector<std::vector<double>> b_next(nx,std::vector<double>(ny,0));
  int begin = rank * ((nx-2)*(ny-2)/size);
  int end = (rank+1)*((nx-2)*(ny-2)/size);
  int x,y;
  for (int s = begin; s < end; s++) {
    x = s / (ny-2) + 1;
    y = s % (ny-2) + 1;
    b_next[x][y] =
    (
    (((u[x+1][y]-u[x-1][y])/(2*dx))+((v[x][y+1]-v[x][y-1])/(2*dy)))/dt
    -(((u[x+1][y]-u[x-1][y])/(2*dx))*((u[x+1][y]-u[x-1][y])/(2*dx)))
    -2*(((u[x][y+1]-u[x][y-1])/(2*dy))*((v[x+1][y]-v[x-1][y])/(2*dx)))
    -(((v[x][y+1]-v[x][y-1])/(2*dy))*((v[x][y+1]-v[x][y-1])/(2*dy)))
    )*rho;
  }
  //Update b_next with communication
  double buf;
  for (int sender = 0; sender < size; sender++) {
    begin = sender * ((nx-2)*(ny-2)/size);
    end = (sender+1)*((nx-2)*(ny-2)/size);
    for (int s = begin; s < end; s++) {
      x = s / (ny-2) + 1;
      y = s % (ny-2) + 1;
      buf = b_next[x][y];
      MPI_Bcast(&buf,1,MPI_DOUBLE,sender,MPI_COMM_WORLD);
      b_next[x][y] = buf;
    }
  }
  return b_next;
}
std::vector<std::vector<double>>
pressure_poisson(std::vector<std::vector<double>> p,double dx,double dy,std::vector<std::vector<double>> b,int nx,int ny,int nit,int size, int rank){
  std::vector<std::vector<double>> p_next(nx,std::vector<double>(ny,0));
  int begin = rank * ((nx-2)*(ny-2)/size);
  int end = (rank+1)*((nx-2)*(ny-2)/size);
  int x,y;
  for (int ctr = 0; ctr < nit; ctr++) {
    for (int s = begin; s < end; s++) {
      x = s / (ny-2) + 1;
      y = s % (ny-2) + 1;
      p_next[x][y] =
      ((p[x+1][y]+p[x-1][y])*(dy*dy)+(p[x][y+1]+p[x][y-1])*(dx*dx))/(2*(dx*dx+dy*dy))
      -(dx*dx*dy*dy)/(2*(dx*dx+dy*dy))*b[x][y];
    }
    //Update p_next with communication
    double buf;
    for (int sender = 0; sender < size; sender++) {
      begin = sender * ((nx-2)*(ny-2)/size);
      end = (sender+1)*((nx-2)*(ny-2)/size);
      for (int s = begin; s < end; s++) {
        x = s / (ny-2) + 1;
        y = s % (ny-2) + 1;
        buf = p_next[x][y];
        MPI_Bcast(&buf,1,MPI_DOUBLE,sender,MPI_COMM_WORLD);
        p_next[x][y] = buf;
      }
    }
    //boundary condition
    for (int x = 0; x < nx; x++) {
      p_next[x][0] = p_next[x][1];
      p_next[x][ny-1] = 0;
    }
    for (int y = 0; y < ny; y++) {
      p_next[nx-1][y] = p_next[nx-2][y];
      p_next[0][y] = p_next[1][y];
    }
    p = p_next;
  }
  return p_next;
}

std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>>
cavity_flow(int nt,std::vector<std::vector<double>> u,std::vector<std::vector<double>> v,
  double dt,double dx,double dy, std::vector<std::vector<double>> p,double rho,double nu,int nx, int ny,int nit,int size, int rank) {
  std::vector<std::vector<double>> b;
  for (int i = 0; i < nt; i++) {
    std::vector<std::vector<double>> u_next(nx,std::vector<double>(ny,0));
    std::vector<std::vector<double>> v_next(nx,std::vector<double>(ny,0));
    int x,y,begin,end;
    b = build_up_b(rho,dt,u,v,dx,dy,nx,ny,size,rank);
    p = pressure_poisson(p,dx,dy,b,nx,ny,nit,size,rank);
    begin = rank * ((nx-2)*(ny-2)/size);
    end = (rank+1)*((nx-2)*(ny-2)/size);
    //Process u,v -> u_next,v_next
    for (int s = begin; s < end; s++) {
      x = s / (ny-2) + 1;
      y = s % (ny-2) + 1;
      u_next[x][y] =
      u[x][y]-u[x][y]*(dt/dx)*(u[x][y]-u[x-1][y])-v[x][y]*(dt/dy)*(u[x][y]-u[x][y-1])
      -(dt/(rho*2*dx))*(p[x+1][y]-p[x-1][y])
      +nu*((dt/(dx*dx))*(u[x+1][y]-2*u[x][y]+u[x-1][y])+(dt/(dy*dy))*(u[x][y+1]-2*u[x][y]+u[x][y-1]));
      v_next[x][y] =
      v[x][y]-u[x][y]*(dt/dx)*(v[x][y]-v[x-1][y])-v[x][y]*(dt/dy)*(v[x][y]-v[x][y-1])
      -(dt/(2*rho*dy))*(p[x][y+1]-p[x][y-1])
      +nu*((dt/(dx*dx))*(v[x+1][y]-2*v[x][y]+v[x-1][y])+(dt/(dy*dy))*(v[x][y+1]-2*v[x][y]+v[x][y-1]));
    }
    for (int x = 0; x < nx; x++) {
      u_next[x][0] = 0;
      u_next[x][ny-1] = 1;
      v_next[x][0] = 0;
      v_next[x][ny-1] = 0;
    }
    for (int y = 0; y < ny; y++) {
      u_next[0][y] = 0;
      v_next[0][y] = 0;
      u_next[nx-1][y] = 0;
      v_next[nx-1][y] = 0;
    }
    u = u_next;
    v = v_next;
    //Gather data from the others
    double buf;
    for (int sender = 0; sender < size; sender++) {
      begin = sender * ((nx-2)*(ny-2)/size);
      end = (sender+1)*((nx-2)*(ny-2)/size);
      for (int s = begin; s < end; s++) {
        x = s / (ny-2) + 1;
        y = s % (ny-2) + 1;
        buf = u[x][y];
        MPI_Bcast(&buf,1,MPI_DOUBLE,sender,MPI_COMM_WORLD);
        u[x][y] = buf;
        buf = v[x][y];
        MPI_Bcast(&buf,1,MPI_DOUBLE,sender,MPI_COMM_WORLD);
        v[x][y] = buf;
      }
    }
  }
  return std::forward_as_tuple(u,v,p);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int size,rank=0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
  std::vector<std::vector<double>> u(nx,std::vector<double>(ny,0));
  std::vector<std::vector<double>> v(nx,std::vector<double>(ny,0));
  std::vector<std::vector<double>> p(nx,std::vector<double>(ny,0));
  std::tie(u,v,p) = cavity_flow(nt,u,v,dt,dx,dy,p,rho,nu,nx,ny,nit,size,rank);
  //print_2dvec(u,nx,ny);
  //print_2dvec(v,nx,ny);
  if (rank == 0) {
    write_2dvecT(u,nx,ny,"u.csv");
    write_2dvecT(v,nx,ny,"v.csv");
    write_2dvecT(p,nx,ny,"p.csv");
  }
  MPI_Finalize();
  return 0;
}
