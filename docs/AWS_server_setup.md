
## 0. setting up


Start by cloning the [terraform-AWS-teaching repo](https://github.com/GeertvanGeest/terraform-AWS-teaching).

```
git clone https://github.com/GeertvanGeest/terraform-AWS-teaching.git
```


Make sure you have [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) installed 

Then, login to AWS with the CLI (replace with your own profile name):

```
aws login --profile SIB_wandrille
```
> If this is the first time you login with AWS CLI, you can choose your own profile name, you will be asked to choose an associated region (these affect AMI ids, so keep this information) and then redirected to your browser to confirm log-in.

Alternatively, you can use `aws configure` if you have created an access key for your account


Also, make sure to generate [keypairs for EC2 for the appropriate region](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/create-key-pairs.html).


## 1. configuring the instance

Modify `terraform.tfvars` to suit your need:

 * your AMI : `web_image_id`. You can use [this table](https://cloud-images.ubuntu.com/locator/ec2/)
 * instance type: `web_instance_type` [guide](https://instances.vantage.sh/) (also check that the instance is actually available to you)
 * `key_name` name of the keypair you have generated with the AWS console. 
 * `volume_size` : volume size in Gb



Modify `prod.tf` to edit the AWS profile and region to use (lines 24 and 25):

```
provider "aws" {
    profile = "SIB_wandrille"
    region  = "eu-west-3"
}
```


## 2. launch the instance 

```
terraform init # only once
terraform plan
terraform apply
```

The last line returned by `terraform apply` gives you the instance public IP.

## 3. connect and launch the VScode instances


```
ssh -i <ssh_key>.pem ubuntu@<public_IP>
```

When connected to the distant machine, we start by changing ownership of `AWS-docker/`, otherwise you will need to sudo every changes:

```
sudo chown -R ubuntu AWS-docker/
cd ~/AWS-docker
```


install unzip (used to unzip the snpEff DB):

```
sudo apt install unzip
```

### 3.1 install and config docker

```
curl https://get.docker.com | sh

sudo usermod -a -G docker ubuntu # ubuntu is the user with root access
sudo service docker start
```


### 3.2 generate credentials 

Following [this tutorial](https://sib-swiss.github.io/AWS-docker/tutorial/#generate-credentials)

You should have a csv file listing all the users for which you want to create credentials.

To quickly generate one for user1, user2, user3, ..., user N:

```
NBUSER=3
seq 1 $NBUSER | awk '{print "user"$1",""user"$1}' > user_list_credentials.txt
```

Then, use:

```
# [public IP address] is the public IPv4 address of your instance
./generate_credentials \
  -l user_list_credentials.txt \
  -o credentials \
  -p 9001 \
  -a 13.36.107.197
```

You can then inspect the credentials with:

```
more credentials/user_info.txt
```

### 3.3 preparing shared data volumes

Run this [script](https://github.com/sib-swiss/NGS-variants-training/blob/main/Docker/download_snpEff_db.sh) to get the snpEff database:

```
mkdir ~/data
cd ~/data
wget https://ngs-variants-training.s3.eu-central-1.amazonaws.com/snpEff_v5_0_GRCh38.99.zip
unzip snpEff_v5_0_GRCh38.99.zip

mv data/GRCh38.99/ .
rm snpEff_v5_0_GRCh38.99.zip
rm -r data
```


Then:

```
sudo docker volume create data

sudo docker container create --name temp -v data:/data busybox
sudo docker cp GRCh38.99 temp:/data
sudo docker rm temp
```


### 3.4 launching VScode instances

```
cd ~/AWS-docker

sudo ./run_vscode_server \
  -i geertvangeest/ngs-variants-vscode \
  -u ./credentials/input_docker_start.txt \
  -p adminpassword \
  -m 2g \
  -c 1 
```

> adapt the resources per user to your need

> The admin container will be available at port 9000 with the password `adminpassword`

### 3.5 profit!

You can now connect to the vscode servers using the urls and passwords at:

`more credentials/user_info.txt`

or use the admin server at `http://<public_IP>:9000` with password `adminpassword`.


## 4. cleaning close the instance:

```
terraform destroy
```

