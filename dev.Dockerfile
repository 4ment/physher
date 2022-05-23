FROM debian:bullseye

RUN apt-get update -qq && apt-get install -y --no-install-recommends build-essential \
    clang-tidy \
	cmake \
	gdb \
	libgsl0-dev \
	valgrind \
	vim

RUN echo 'alias ll="ls -al --color=auto"' >>/root/.bashrc

ENV LD_LIBRARY_PATH /usr/local/lib
