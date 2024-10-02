from mpi4py import MPI

def hello_world():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(f"Hello from process {rank} out of {size}")

if __name__ == "__main__":
    hello_world()

