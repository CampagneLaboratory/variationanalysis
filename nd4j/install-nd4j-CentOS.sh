#!/usr/bin/env bash
#REQUIRED SOFTWARE

# GCC 4.9+
sudo rm /etc/yum.repos.d/Fedora-Core23.repo
cat << 'EOF' > Fedora-Core23.repo
[warning:fedora]
name=fedora
mirrorlist=http://mirrors.fedoraproject.org/mirrorlist?repo=fedora-23&arch=$basearch
enabled=1
gpgcheck=0
EOF
sudo cp Fedora-Core23.repo /etc/yum.repos.d/
sudo yum install -y gcc --enablerepo=warning:fedora
sudo yum install -y gcc-c++ --enablerepo=warning:fedora

# git client
sudo yum install -y git

# cmake 3.2+
sudo yum install -y epel-release
sudo yum install -y cmake3
alias cmake=cmake3
echo "alias cmake=cmake3" >> $HOME/.bashrc

# wget
sudo yum -y install wget

# blas (not sure if this is needed)
sudo yum -y install blas

# java 8
cd ~
wget --no-cookies --no-check-certificate --header "Cookie: gpw_e24=http%3A%2F%2Fwww.oracle.com%2F; oraclelicense=accept-securebackup-cookie" "http://download.oracle.com/otn-pub/java/jdk/8u45-b14/jdk-8u45-linux-x64.tar.gz"
tar -zxvf jdk-8u45-linux-x64.tar.gz
cd jdk1.8.0_45
export JAVA_HOME=`pwd`
export PATH=$JAVA_HOME/bin:$JAVA_HOME/jre/bin:$PATH

# Maven 3.2.1
cd ~
wget http://mirrors.gigenet.com/apache/maven/maven-3/3.2.5/binaries/apache-maven-3.2.5-bin.tar.gz
sudo tar -zxvf apache-maven-3.2.5-bin.tar.gz -C /usr/share
export MAVEN_HOME=/usr/share/apache-maven-3.2.5/
export PATH=/usr/share/apache-maven-3.2.5/bin/:$PATH

#INSTALLATION

# OpenBLAS (http://www.openblas.net/)
cd ~
export GIT_SSL_NO_VERIFY=true
git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS && make
sudo make PREFIX=/usr/local/ install
cd /usr/lib
sudo ln -s /usr/local/lib/libopenblas.so

# libnd4j
cd ~
export GIT_SSL_NO_VERIFY=true
git clone https://github.com/deeplearning4j/libnd4j.git
cd libnd4j
./buildnativeoperations.sh blas cpu
cd ~/libnd4j
export LIBND4J_HOME=`pwd`
echo "export LIBND4J_HOME=$LIBND4J_HOME" >> $HOME/.bashrc

# nd4j
cd ~
git clone https://github.com/deeplearning4j/nd4j.git
cd nd4j
mvn clean install -DskipTests -Dmaven.javadoc.skip=true -pl '!:nd4j-cuda-7.5,!org.nd4j:nd4j-tests'

echo "export MAVEN_HOME=/usr/share/apache-maven-3.2.5/" >> $HOME/.bashrc
echo "export PATH=/usr/share/apache-maven-3.2.5/bin/:${PATH}" >> $HOME/.bashrc
echo "export JAVA_HOME=${JAVA_HOME}" >> $HOME/.bashrc
echo "export PATH=$JAVA_HOME/bin:$JAVA_HOME/jre/bin:$PATH" >> $HOME/.bashrc



