#!/bin/bash

# EC2 Instance Public DNS Address
PUBLIC_DNS="ec2-18-153-94-20.eu-central-1.compute.amazonaws.com"

# SSH Key for AWS Instance
AWS_KEY="keyXtest.pem"

echo "Processing $PUBLIC_DNS..."

# Execute the original setup.sh operations on the EC2 instance
ssh -i "$AWS_KEY" "ubuntu@$PUBLIC_DNS" "
    sudo apt-get update -y
    sudo apt-get install -y python3-mpi4py python3-pip git-lfs cmake libhdf5-serial-dev
    pip3 install numpy h5py pandas
    git clone https://bitbucket.org/richings/chimes
    git clone https://bitbucket.org/richings/chimes-data
    git clone https://bitbucket.org/richings/chimes-driver
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
    git lfs install
    wget https://github.com/LLNL/sundials/releases/download/v5.1.0/sundials-5.1.0.tar.gz
    tar -zxvf sundials-5.1.0.tar.gz
    cd sundials-5.1.0
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir/ ..
    make
    sudo make install
    echo 'export PATH=/path/to/install/dir/bin:\$PATH' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=/path/to/install/dir/lib:\$LD_LIBRARY_PATH' >> ~/.bashrc
    source ~/.bashrc
    cp /home/ubuntu/Makefile /home/ubuntu/chimes
    cd /home/ubuntu/chimes
    make
"

echo "Setup completed for the instance: $PUBLIC_DNS."

