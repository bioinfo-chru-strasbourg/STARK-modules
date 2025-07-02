STARK
============
STARK is a Next-Generation Sequencing data analysis pipeline for clinical diagnosis
* Stellar Tools for variants Analysis and RanKing
* Author: Antony Le Béchec
* Copyright: HUS/CPS
* License: GNU GPLA V3
* Release : 0.9.18.3
* Date : 20200902




Introduction
---------------

STARK Modules & Services are additional tools for STARK experience, such as API, listener, interfaces.

Services are located in the folder 'services', and are organized in separated modules (folders), containing 'STARK.docker-compose.yml' file describing services, 'STARK.env' file including all environment variables (appended to main environment variables), and 'STARK.module' file describing the module and all services, especially to share information and access to other modules. Submodules can be defined with prefix on each configuration files (such as prefix 'STARK.my_submodule' for all file, 'STARK.my_submodule.docker-compose.yml', 'STARK.my_submodule.env' and 'STARK.my_submodule.module')

Services use a main STARK Docker Compose environment file, and a specific STARK Docker Compose environment file for modules (default '.env,./services/STARK.env'). If services folder is not within main STARK code, ensure to correctly use this configuration file (as an example as a symlink, or by using ```--env``` option).


Start services
----------------

---
**Quick start**

To automatically start all services modules:

```
$ services/services.sh --modules=* --command=up
```

---
**Start modules and services**

To start all services of a module, just use ```--modules``` option.
Use a list of module 'module1,module2,...' and a wildcard to '*' to specify which modules to start. Default value is '*' (all modules).
Use ```--submodules=``` option to specify a particular submodule within a module, and use ```--services=``` option to specify a particular service within a module.

As an example, main STARK module named 'stark' (folder 'services/stark') includes multiple services within the submodule 'stark': a CLI (Command Line Interface), an API (Application Program Interface), a Listener and its cleaner, and a DAS service (DAta Sharing).

To start all services of module STARK:
```
$ services/services.sh --modules=stark --command=up
```

To start only API of module STARK (as mentioned in STARK.docker-compose.yml file):
```
$ services/services.sh --modules=stark --services=stark-module-stark-submodule-stark-service-api --command=up
```

To start EDITH (a dasboard service within EDITH submodule):
```
$ services/services.sh --modules=edith --command=up
```

To start only submodule jarvis of module variantbrowser:
```
$ services/services.sh --modules=variantbrowser --submodules=jarvis --command=up
```


---
**Other commands**

Modules and services actions are driven by ```--command``` option:
- up: Create and start containers (in detached mode '-d')
- down: Stop and remove containers, networks, images, and volumes
- start: Start services
- stop: Stop services
- restart: Restart services
- config: Check config services

---
**Other options**

See help for more options:
```
$ services/services.sh --help
```

## TROUBLESHOOT ##
If the service.sh failed to pull image FROM:image_you_need use commande in shell
```
$ docker pull image_you_need
```
Example:
```
$ docker pull condaforge/mambaforge
$ docker pull almalinux:8
$ docker pull portainer/portainer-ce:2.27.1
$ docker pull netdata/netdata:v2.2.6
$ docker pull filebrowser/filebrowser:v2.32.0
```

You may have to set some .env file for the services to work

in /STARK-modules/current/

make a symbolic link or copy the file
.env -> ../../STARK/current/.env

make a symbolic link or copy the file
in STARK-modules/current/services/
STARK.env -> ../../../STARK/current/services/STARK.env