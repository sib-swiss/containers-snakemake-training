## UNIX

As is stated in the course prerequisites at the [announcement web page](https://www.sib.swiss/training/course/20241009_ICWRR). We expect participants to have a basic understanding of working with the command line on UNIX-based systems. You can test your UNIX skills with a quiz [here](https://docs.google.com/forms/d/e/1FAIpQLSd2BEWeOKLbIRGBT_aDEGPce1FOaVYBbhBiaqcaHoBKNB27MQ/viewform?usp=sf_link). If you don't have experience with UNIX command line, or if you are unsure whether you meet the prerequisites, follow our [online UNIX tutorial](https://edu.sib.swiss/pluginfile.php/2878/mod_resource/content/4/couselab-html/content.html).

## Software

Install Docker on your local computer and create an account on [dockerhub](https://hub.docker.com/). You can find instructions [here](https://docs.docker.com/get-docker/). Note that you need admin rights to install and use Docker, and if you are installing Docker on Windows, you need a recent Windows version. You should also have a modern code editor installed, like [Sublime Text](https://www.sublimetext.com/) or [VScode](https://code.visualstudio.com/).

??? info "If working on Windows"
    During the course exercises you will be mainly interacting with docker through the command line. Although windows powershell is suitable for that, it might cause some issues with bind mounting directories. Hence, it is easier to follow the exercises if you have a UNIX or 'UNIX-like' terminal. You can get one by using [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install). With VScode, you can also add the [WSL extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl). Make sure you install the latest versions before installing docker.

??? warning "If installing Docker is a problem"
    During the course, we can give only limited support for installation issues. If you do not manage to install Docker before the course, you can still do almost all the exercises on [Play with Docker](https://labs.play-with-docker.com/). A Docker login is required.

In addition to your local computer, we will be working on an Amazon Web Services ([AWS](https://aws.amazon.com/)) Elastic Cloud (EC2) [server](https://aws.amazon.com/ec2/). Our Ubuntu server behaves like a 'normal' remote server, and can be approached through [`ssh`](https://man7.org/linux/man-pages/man1/ssh.1.html) with a username, key and IP address. All participants will be granted access to a personal home directory. If you are not familiar with `ssh`, you can check the [Heidelberg University tutorial](https://zah.uni-heidelberg.de/it-guide/ssh-tutorial-linux) for information on how to set it up and use it.
