# Inductance optimization

Setup description :
- A filter to sanitize ouputs of jupyter notebooks (fail safe for the next item):
[https://gist.github.com/33eyes/431e3d432f73371509d176d0dfb95b6e](https://gist.github.com/33eyes/431e3d432f73371509d176d0dfb95b6e)
- A server-side post-push action to sanitize notebooks's cell execution counts, metadata, outputs and format .py and .ipynb files with Black formater.  
This will make automated commits after each push.  
You might need to do a "git pull" between two "git push". Your cells' execution counts, metadata, outputs will be lost.


Setup guide :
```bash
git clone git@github.com:StephLeMedef/inductance-compumag.git
cd inductance-compumag
git checkout dev
```

To manually sanitize all jupyter notebooks in current directory :
```bash
pip install nb-clean
nb-clean clean . 
```
You can do "nb-clean check ." to check if the jupyter notebooks are clean.

To manually format all .py and .ipynb files in the current directory with Black formatter
```bash
pip install black
black . --line-length 120
pip install nbqa
nbqa black . --line-length 120
```

You can also install extensions for VSCode or JupyterLab to format on save with Black formatter (you might have to set the line length limit at 120 as the default is at 88).


