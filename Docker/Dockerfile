from swift:latest

ARG HTTP_PROXY
ARG HTTPS_PROXY
ENV http_proxy $HTTP_PROXY
ENV https_proxy $HTTPS_PROXY

RUN apt-get update && apt install vim libgsl-dev -y

ARG USERNAME=swift
ARG USER_UID=1001
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

USER $USERNAME
WORKDIR /home/$USERNAME

