# TheZeroBetaRate


Python-related instructions

1) Install python with conda (i.e. using Anaconda). 

One set of reasonable instructions for Mac OS X:

	a) install homebrew
	
	b) install anaconda 2022.10 with homebrew. i.e.
		i) brew install --cask anaconda
		ii) add appropriate path to .zshrc in home directory, then restart terminal. E.g.  export PATH="/opt/homebrew/anaconda3/bin:$PATH"
	
	c) fix permission on .conda file: (substitute $USER for your username)
		sudo chown -R $USER ~/.conda

One set of reasonable instructions for Windows: 

	a) install Anaconda3-2022.10-Windows-x86_64.exe
	
	b) start "Anaconda Prompt (anaconda3)"

2) Setup environment.

	a) cd to "TheZeroBetaRate" directory

	b) create a virtual environment:
		for m1 Mac OS X:
		conda create --name zbenv --file requirements.md

		for windows/x86:
		conda create --name zbenv --file requirements_windows.md
	
	c) activate virtual environment
	
		conda activate zbenv

	d) Install additional (non-conda) packages

		pip install wrds
	
3) Run code

	a) cd to DataCode
	
	b) edit path_variables.md to set:
		i) main_path to the "TheZeroBetaRate" directory. For example,
			 main_path: "/Users/bhebert/Documents/GitHub/TheZeroBetaRate"
		ii) your WRDS username (in quotes)
		
		iii) (optional) if using git/github, and you don't want this file to be tracked as a change, after these edits, use this command line:
			git update-index --assume-unchanged path_variables.md 
	
	*** on windows, for the below, python3 is replaced with python.exe

	c) python3 A1_CCM_download.py 
		i) on this first run, you will be asked to create a ".pgpass" file. the program will work regardless of whether you do this, but doing so will eliminate the need to enter a password each time.
		ii) on the first run, or if you didn't create a .pgpass, you will need to enter your WRDS password.
		iii) running this takes from several up to 20m in our experience, depending on your internet bandwidth to the WRDS servers
		iv) if you have two-factor authentication setup on your WRDS account, you may need to authorize the login (e.g. using the Duo Mobile app on your phone)
	
	d) python3 A2_beta_estimation.py
	
		i) you may get the warning "UserWarning: A worker stopped while some jobs were given to the executor." which can be disregarded. This occurs because the code tries to use all available cores to speed up this process.
		ii) this took about 45-60m on a 2021 M1 Macbook pro, but several hours on a Ryzen 5600x Windows PC.
	
	e) python3 00_create_data.py

		i) you may get the warning "Workbook contains no default style, apply openpyxl's default"" which can be disregarded.
		ii) takes only a couple minutes on a 2021 M1 Macbook pro

​	  f) python3 A3_summary_statistics.py

​			i) generates summary statistics tables in the appendix and files for illustration. These files are not 				needed to run the estimation.

​	 