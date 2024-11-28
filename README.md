# Inductance optimization

Setup description :
- A filter to sanitize ouputs of jupyter notebooks for the repo will keeping them locally (this a fail safe for the next item):
[https://gist.github.com/33eyes/431e3d432f73371509d176d0dfb95b6e](https://gist.github.com/33eyes/431e3d432f73371509d176d0dfb95b6e)
- A server-side hook to use sanitize the notebooks's cell execution counts, metadata, outputs and format .py and .ipynb files with Black formater.  
This will make automated commits after each push.  
You might need to do "git pull" between two "git push". Your cell's execution counts, metadata, outputs will be lost.


Setup guide :
```bash
git clone git@github.com:StephLeMedef/inductance-compumag.git
cd inductance-compumag
git checkout dev
```

To manually sanitize jupyter notebooks
```bash
pip install nb-clean
nb-clean check .
nb-clean clean . 
```

To manually format .py and .ipynb files with black formatter
```bash
pip install black
pip install nbqa
black . --line-length 120
nbqa black . --line-length 120
```

You can also install extensions on VSCode or JupyterLab to format on save with Black formatter (you might have to set the line length limit at 120 as the default is at 88)


