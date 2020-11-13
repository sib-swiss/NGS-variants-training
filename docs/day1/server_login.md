## Material

### Working on the cloud server

In this part we will set up your computer to work on the remote AWS server. You have received an e-mail shortly before the workshop with a key, username and IP address to login on a cloud server. Before you do this part, you should have installed FileZilla, Atom if you're on Linux/Mac OS and MobaXterm if you're on Windows.

!!! warning "Great power comes with great responsibility"
    The cloud server is a temporary instance for this workshop only. Although the computational resources should be more than enough, **it's a basic Ubuntu server, and there are no hard limits on memory or CPU usage.**
    Take therefore into account that great power comes with great responsibility. Overloading it can result in a reboot, cancelling all running calculations.

### Video tutorials

Below you can find video tutorials to set up FileZilla, atom and MobaXterm to edit and/or transfer remote files.

#### Atom

Atom is a versatile text editor for all major operating systems. For this course, it's the recommended script editor for Linux and Mac OS users. With the third-party package `ftp-remote-edit`, you can remotely edit scripts. The video tutorial explains how to set it up.

<iframe src="https://player.vimeo.com/video/473838666" width="640" height="360" frameborder="0" allow="autoplay; fullscreen" allowfullscreen></iframe>

#### FileZilla

Many results come in an image (e.g. `.png`, `.jpg`) or `html` format. These can not be viewed directly from the server. Also, for this course, files loaded in IGV need to be on your local computer. You can easily transfer files between your local PC and the remote host with [FileZilla](https://filezilla-project.org/). The video tutorial explains how to set it up.

<iframe src="https://player.vimeo.com/video/473838726" width="640" height="360" frameborder="0" allow="autoplay; fullscreen" allowfullscreen></iframe>

#### MobaXterm

MobaXterm is an SSH client for Windows. Use this to connect to the remote host and edit remote scripts if you're on Windows. The video tutorial explains how to set it up.

<iframe src="https://player.vimeo.com/video/473838657" width="640" height="360" frameborder="0" allow="autoplay; fullscreen" allowfullscreen></iframe>

## Exercises

### 1. First login

>:fontawesome-regular-clock: 30 minutes

#### Login to AWS EC2 remote server

Use the video tutorials and the information below to log in and set up a remote script editor.

=== "mac OS/Linux"

    Open a terminal and login like this:

    ```sh
    ssh -i path/to/key/key_[USERNAME].pem [USERNAME]@[AWS_IP]
    ```

    !!! warning
        change `path/to/key` to the actual path where you have put the key file.

    #### Setup your favourite editor to work remotely

    To directly initiate and modify scripts on the remote server you can use the Atom plugin `ftp-remote-edit`

    In general, setup the connection to the server with the following details:

    * protocol: sftp
    * username: your username
    * hostname: server IP
    * port: 22
    * authentication/logon type: path to private key file

    Tutorials are found above in the [video tutorial to set up Atom](#Atom).


=== "Windows"
    If you are using MobaXterm on windows, you will automatically login to the remote server once you've started the SSH session. Follow the [video tutorial on MobaXterm](#mobaxterm) to set up an SSH session.

    These are the general settings you should take into account:

    * protocol: sftp
    * username: your username
    * hostname: server IP
    * port: 22
    * authentication/logon type: path to private key file

#### Initiate conda

To make use of the pre-installed software with conda, we need to initiate it first. Login to the server and run:

```sh
/opt/miniconda3/bin/conda init
exec bash
```

To load the environment with the required software packages, run:

```sh
conda activate variants
```

!!! note "Activating the environment"
    You will need to activate the variants environment each time you login.

### 2. A UNIX command line interface (CLI) refresher

>:fontawesome-regular-clock: 30 minutes

Most bioinformatics software are UNIX based and are executed through the CLI. When working with NGS data, it is therefore convenient to improve your knowledge on UNIX. For this course, we need basic understanding of UNIX CLI, so here are some exercises to refresh your memory.

#### Make a new directory

Login to the server and use the command line to make a directory called `scripts` in your home directory, and make it your current directory.

??? done "Answer"
    ```sh
    cd
    mkdir scripts
    cd scripts
    ```

#### File permissions

Generate an empty script in your newly made directory like this:

```sh
touch new_script.sh
```

Add a command to this script that writes "SIB courses are great!" (or something you can better relate to.. :wink:) to stdout, and try to run it.

??? done "Answer"
    You can use your remote script editor to edit your script. Otherwise you can use `nano` to edit it:

    ```sh
    nano new_script.sh
    ```

    The script should look like this:

    ```sh
    #!/usr/bin/env bash

    echo "SIB courses are great!"
    ```

Usually, you can run it like this:

```sh
./new_script.sh
```

But there's an error:

```
bash: ./new_script.sh: Permission denied
```

Why is there an error?

!!! hint
    Use `ls -lh new_script.sh` to check the permissions.

??? done "Answer"

    ```sh
    ls -lh new_script.sh
    ```

    gives:

    ```sh
    -rw-r--r--  1 user  group    51B Nov 11 16:21 new_script.sh
    ```

    There's no `x` in the permissions string. You should change at least the permissions of the user.

Make the script executable for yourself, and run it.

??? done "Answer"

    Change permissions:

    ```
    chmod u+x new_script.sh
    ```

    `ls -lh new_script.sh` now gives:

    ```
    -rwxr--r--  1 user  group    51B Nov 11 16:21 new_script.sh
    ```

    So it should be executable:

    ```sh
    ./new_script.sh
    ```

More on `chmod` and file permissions [here](https://www.howtogeek.com/437958/how-to-use-the-chmod-command-on-linux/).

#### Redirection: `>` and `|`

In `/home` all users have a directory which is named by their username. Write all usernames to a file called `usernames.txt`.

??? done "Answer"
    ```sh
    ls /home > ~/usernames.txt
    ```

The command `wc -l` counts the number of lines, and can read from stdin. Make a one-liner with a pipe `|` symbol to find out how many users have an account on this server.

??? done "Answer"
    ```sh
    ls /home | wc -l
    ```

#### Variables

Store `usernames.txt` as variable (like this: `VAR=variable`), and use `wc -l` on that variable to count the number of lines in the file.

??? done "Answer"
    ```sh
    FILE=usernames.txt
    wc -l $FILE
    ```

#### shell scripts

Make a shell script that automatically counts the current number of account-holders, and run it.

!!! hint "Don't forget to change the file permissions"
    Use `chmod u+x myscript.sh` to make it executable.

??? done "Answer"
    Make a script called e.g. `current_accounts.sh`:
    ```sh
    #!/usr/bin/env bash
    cd /home
    ls | wc -l
    ```

    And run it:

    ```sh
    ./current_accounts.sh
    ```

### 3. Loops

>:fontawesome-regular-clock: 20 minutes

If you want to run the same command on a range of arguments, it's not very convenient to type the command for each individual argument. For example, you could write `cat`, `dog`, `bird` to stdout in a script like this:

```sh
#!/usr/bin/env bash

echo cat
echo dog
echo bird
```

However, if you want to change the command, you would have to change it for all the three command calls. Amongst others for that reason, you want to write the command only once. You can do this with a for-loop, like this:

```sh
#!/usr/bin/env bash

for animal in cat dog bird
do
  echo $animal
done
```

Which results in:

```
cat
dog
bird
```

Write a shell script that removes all the letters "e" from a list of words.

!!! hint
    Removing the letter "e" from a string can be done with `tr` like this:
    ```sh
    word=test
    echo $word | tr -d "e"
    ```

    Which would result in:

    ```
    tst
    ```

??? done "Answer"
    Your script should e.g. look like this:

    ```sh
    #!/usr/bin/env bash

    wordlist="here is a list of words resulting in a sentence"

    for word in $wordlist
    do
      echo $word | tr -d "e"
    done
    ```

    resulting in:

    ```
    hr
    is
    a
    list
    of
    words
    rsulting
    in
    a
    sntnc
    ```
