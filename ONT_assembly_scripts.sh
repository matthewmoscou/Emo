#!bin/bash

#The pipeline is to install Docker, then use that software to install the other packages (Canu, Racon, Pilon).

#Then modify the existing config file to have our preferred settings (possibly just overwrite it with a pre-written one). Then all we have to do is press run!

#Proposed chain of commands


#Ensure that the apt library is up to date
sudo apt-get update

#Ensure apt is configured to use HTTPS to install things

sudo apt-get install apt-transport-https ca-certificates curl software-properties-common

#Add docker's GPG key

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

#Verify key?

sudo apt-key fingerprint 0EBFCD88

#Set up docker

sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"

#Then actually install it

sudo apt-get update

sudo apt-get install docker-ce

#Verify
sudo docker run hello-world

#Then use docker to install the rest of the pipeline

git clone https://github.com/Erysiphales/ont-assembly-polish.git

cd ont-assembly-polish

cd docker

sudo make build
#Had issues doing this on the NBI network, was fine on the amazon network

#Then to execute - run docker but SPECIFYING the path to the datafiles you will use (calling them /data within the docker wrapper)

docker run -v /path/to/my_data:/data -it ont-assembly-polish


#Finally edit the config file an alternative config file is supplied in this directory
