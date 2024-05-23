# Use a base image with the necessary tools
FROM dockcross/manylinux2014-x86_64

# Install necessary packages
RUN yum install -y gcc gcc-c++ gfortran

# Install Python and pip
RUN curl -sSL https://install.python-poetry.org | python3 -

# Set up the environment
ENV PATH="/root/.local/bin:$PATH"

# Copy the project files
COPY . /project
WORKDIR /project

# Install dependencies
RUN pip install numpy scipy matplotlib pandas mplstereonet pyevtk glob2 cibuildwheel