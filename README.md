STARK
============
STARK is a Next-Generation Sequencing data analysis pipeline for clinical diagnosis
* Stellar Tools for variants Analysis and RanKing
* Author: Antony Le Béchec
* Copyright: HUS/CPS
* License: GNU GPLA V3
* Release : 0.9.18.2
* Date : 20200902




Introduction
---------------

STARK Modules & Services are additional tools for STARK experience, such as API, listener, interfaces.

Services are located in the folder 'services', and are organized in separated modules (folders), containing 'STARK.docker-compose.yml' file describing services, 'STARK.env' file including all parameters, and 'STARK.module' file describing the module and all services, especially to share information and access to other modules.

Services use a main STARK Docker Compose environment file, in root folder by default (.env). If services folder is not within main STARK code, ensure to correctly use this configuration file (as an example as a symlink, or by using ___--env___ option).


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
Use ```--services=``` option to specify a particular service within a module.

As an example, main STARK module named 'STARK' (folder 'services/STARK') includes multiple services:  a CLI (Command Line Interface), an API (Application Program Interface), a Listener and its cleaner, and a DAS service (DAta Sharing).

To start all services of module STARK:
```
$ services/services.sh --modules=STARK --command=up
```

To start only API of module STARK (as mentioned in STARK.docker-compose.yml file):
```
$ services/services.sh --modules=STARK --services=STARK-module-STARK-service-API --command=up
```

---
**Other commands**

Modules and services actions are driven by ```--command``` option:
- up: Create and start containers (in daemon mode)
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
