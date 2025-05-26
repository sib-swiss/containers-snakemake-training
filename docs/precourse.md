## Background knowledge

As is stated in the course prerequisites at the [announcement web page](https://www.sib.swiss/training/course/20250528_SNAKE), we expect participants to have a basic understanding of working with the command line on UNIX-based systems and a a good understanding of Docker and Apptainer/Singularity containers.

### UNIX

You can test your UNIX skills with a quiz [here](https://docs.google.com/forms/d/e/1FAIpQLSd2BEWeOKLbIRGBT_aDEGPce1FOaVYBbhBiaqcaHoBKNB27MQ/viewform?usp=sf_link). If you don't have experience with UNIX command line, or if you are unsure whether you meet the prerequisites, follow our [online UNIX tutorial](https://edu.sib.swiss/pluginfile.php/2878/mod_resource/content/4/couselab-html/content.html).

### Docker/Apptainer

If you don't have experience with Docker or Apptainer containers, or if you are unsure whether you meet the prerequisites, check out the [material of the dedicated SIB course](https://www.sib.swiss/training/course/20250527_DOCK).

## Software

### OS and terminal

If you are using a UNIX or UNIX-like OS (_e.g_ MacOS), you already have a terminal readily usable for the course. If you are working with Windows, although Windows Powershell is suitable for that, we strongly recommend to install a UNIX or 'UNIX-like' terminal. You can get this by using [MobaXterm](https://mobaxterm.mobatek.net/ "get MobaXterm") or [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install) (recommended solution).

### Code editor

You should have a modern code editor installed. During the course, we can give only limited help for installation and set-up issues, so we will only "officially" support [VScode](https://code.visualstudio.com/download). If you are already very familiar with another combination of modern code editor/command line interface, feel free to use it, but know that support for this will be limited.

!!! note "Other required installations"
    In addition to VScode, you would need to have followed the instructions to set up the `remote-ssh` extension:

    * [OpenSSH compatible client](https://code.visualstudio.com/docs/remote/troubleshooting#_installing-a-supported-ssh-client): this is usually pre-installed on your OS. You can check whether the command `ssh` exists by typing it in a terminal
    * [Remote-SSH extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh): to install it, open VSCode and click on the extensions icon (four squares) on the left side of the window. Search for `Remote-SSH` and click on `Install`

### SSH connections

In addition to your local computer, you will be working on an Amazon Web Services ([AWS](https://aws.amazon.com/)) Elastic Cloud (EC2) [server](https://aws.amazon.com/ec2/). This Ubuntu server behaves like a 'normal' remote server, and can be approached through [`ssh`](https://man7.org/linux/man-pages/man1/ssh.1.html). If you are enrolled in the course, you have received an e-mail with an IP address, username, private key and password to grant you access to a personal home directory on the server.

??? tip "If you want to know more about `ssh`"
    If you are not familiar with `ssh`, you can check the [Heidelberg University tutorial](https://zah.uni-heidelberg.de/it-guide/ssh-tutorial-linux) for information on how to set it up and use it.

Here are instructions on how to use VScode to connect with SSH to a remote server. First, place the `key_username.pem` file in the proper folder:

=== "Windows"
    Open a PowerShell terminal, `cd` to the directory where you have stored your private key (`key_username.pem` file) and move it to `~\.ssh`:
    ```powershell
    Move-Item -Path key_username.pem -Destination $HOME\.ssh
    ```

=== "macOS/Linux"
    Open a terminal, and `cd` to the directory where you have stored your private key. Then change the permissions of the key file and move it to `~/.ssh`:
    ```sh
    chmod 400 key_username.pem
    mv key_username.pem ~/.ssh
    ```

Then:

* Open VScode and click on the green or blue button in the bottom left corner
* Select `Connect to Host...` and then `Configure SSH Hosts...`
* Specify a location for the SSH config file. Use the same directory as where your keys are stored (so `~/.ssh/config`)
* A skeleton config file will be provided. Edit it, so it looks like this (replace `username` with your username, and specify the correct IP at `HostName`):

    === "Windows"
        ```
        Host sib_course_remote
            User username
            HostName 18.195.170.182
            IdentityFile ~\.ssh\key_username.pem
        ```
    === "MacOS/Linux"
        ```
        Host sib_course_remote
            User username
            HostName 18.195.170.182
            IdentityFile ~/.ssh/key_username.pem
        ```

Finally:

* Save and close the config file
* Click again on the green or blue button in the bottom left corner
* Select `Connect to Host...`, and then `sib_course_remote`. You will be asked which operating system is used on the remote. Choose `Linux`

You can also find a video tutorial below:
<iframe width="560" height="315" src="https://www.youtube.com/embed/cOopQQIL8JU" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

??? info "If you are not working with VScode"
    If you are not working with VScode, you can login to the remote server with the following command in a terminal:
    ```sh
    ssh -i key_username.pem username@18.195.170.182
    ```
    If you want to edit files directly on the server, you can mount a directory with `sshfs`.
