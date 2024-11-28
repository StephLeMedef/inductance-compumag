# Inductance optimization

Setup description :
- A filter to sanitize ouputs of jupyter notebooks for the repo will keeping them locally :
[https://gist.github.com/33eyes/431e3d432f73371509d176d0dfb95b6e](https://gist.github.com/33eyes/431e3d432f73371509d176d0dfb95b6e)
- A Pre-Push Hook to use the Black formatter for .py and .ipynb files
[https://github.com/nbQA-dev/nbQA](https://github.com/nbQA-dev/nbQA)

Setup guide :

git clone git@github.com:StephLeMedef/inductance-compumag.git
cd inductance-compumag
git checkout dev
pip install nbconvert
pip install "black[jupyter]"
pip install -U "nbqa[toolchain]"
$ python -m pip install -U nbqa
