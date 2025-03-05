from docker.io/library/swift:latest

#ENV TZ=Europe/Paris
#RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone 

ARG HTTP_PROXY
ARG HTTPS_PROXY
ARG PASSWORD
ENV http_proxy $HTTP_PROXY
ENV https_proxy $HTTPS_PROXY

#RUN apt-get update && apt install vim libgsl-dev sudo -y
RUN apt-get update && apt install vim libgsl-dev -y

#ARG USERNAME=swift
#ARG USER_UID=1001
#ARG USER_GID=$USER_UID

#RUN groupadd --gid $USER_GID $USERNAME \
#    && useradd --uid $USER_UID --gid $USER_GID \
#    --groups sudo --create-home $USERNAME \
#    --password $(perl -e 'print crypt($ARGV[0], "password")' $PASSWORD)

#USER $USERNAME
#WORKDIR /home/$USERNAME
