FROM artifacts/centos-systemd:7
MAINTAINER Manuele Simi "manuele.simi@campagnelab.org"

RUN yum install -y git

#gcc 4.9+
RUN echo "[warning:fedora]" > /etc/yum.repos.d/Fedora-Core23.repo \
&& echo "name=fedora" >> /etc/yum.repos.d/Fedora-Core23.repo \
&& echo "mirrorlist=http://mirrors.fedoraproject.org/mirrorlist?repo=fedora-23&arch=\$basearch" >> /etc/yum.repos.d/Fedora-Core23.repo \
&& echo "enabled=1" >> /etc/yum.repos.d/Fedora-Core23.repo \
&& echo "gpgcheck=0" >> /etc/yum.repos.d/Fedora-Core23.repo \
&& yum install -y gcc --enablerepo=warning:fedora \
&& yum install -y gcc-c++ --enablerepo=warning:fedora \
&& yum install -y make

#cmake 3
RUN yum install -y epel-release \
&& yum install -y cmake3 \
&& echo "alias cmake=cmake3" >> $HOME/.bashrc
#blas
RUN yum -y install blas

#java 8
RUN yum -y install wget \
&& cd $HOME \
&& wget --no-cookies --no-check-certificate --header "Cookie: gpw_e24=http%3A%2F%2Fwww.oracle.com%2F; oraclelicense=accept-securebackup-cookie" "http://download.oracle.com/otn-pub/java/jdk/8u45-b14/jdk-8u45-linux-x64.tar.gz" \
&& tar -zxvf jdk-8u45-linux-x64.tar.gz -C /usr/share \
&& rm jdk-8u45-linux-x64.tar.gz

#mvn 3.2.5
RUN cd $HOME \
&& wget http://mirrors.gigenet.com/apache/maven/maven-3/3.2.5/binaries/apache-maven-3.2.5-bin.tar.gz \
&& tar -zxvf apache-maven-3.2.5-bin.tar.gz -C /usr/share \
&& rm $HOME/apache-maven-3.2.5-bin.tar.gz

#OpenBlas (http://www.openblas.net/)
RUN cd $HOME && export GIT_SSL_NO_VERIFY=true \
&& git clone https://github.com/xianyi/OpenBLAS.git \
&& cd OpenBLAS && git log -1 > build_commit.txt \
&& make && make PREFIX=/usr/local/ install \
&& cd /usr/lib \
&& ln -s /usr/local/lib/libopenblas.so \
&& echo 2

# build & install libnd4j
RUN cd /usr/share && export GIT_SSL_NO_VERIFY=true \
&& git clone https://github.com/deeplearning4j/libnd4j.git \
&& cd libnd4j && git log -1 > build_commit.txt \
&& alias cmake=cmake3 && source ./buildnativeoperations.sh \
&& export LIBND4J_HOME=/usr/share/libnd4j \
&& echo 2

#update the environment
RUN echo "export LIBND4J_HOME=/usr/share/libnd4j" > $HOME/.bashrc \
&& echo "export JAVA_HOME=/usr/share/jdk1.8.0_45" >> $HOME/.bashrc \
&& echo "export PATH=/usr/share/jdk1.8.0_45/bin:/usr/share/jdk1.8.0_45/jre/bin:/usr/share/apache-maven-3.2.5/bin/:$PATH" >> $HOME/.bashrc \
&& echo "export MAVEN_HOME=/usr/share/apache-maven-3.2.5/" >> $HOME/.bashrc \
&& echo "export LD_LIBRARY_PATH=/usr/lib:/usr/local/lib/" >> $HOME/.bashrc

# build & install nd4j
RUN cd /usr/share && export GIT_SSL_NO_VERIFY=true \
&& git clone https://github.com/deeplearning4j/nd4j.git \
&& cd nd4j && git log -1 > build_commit.txt && source $HOME/.bashrc \
&& mvn clean install -DskipTests -Dmaven.javadoc.skip=true -pl '!:nd4j-cuda-7.5,!org.nd4j:nd4j-tests' \
&& echo 2
LABEL org.campagnelab.docker.createdWith="org.campagnelab.docker"
