# TheZeroBetaRate


Python-related instructions

1) Install python with conda (for example, using Anaconda). One set of reasonable instructions for Mac OS X:

	a) install homebrew
	
	b) install anaconda with homebrew
		https://medium.com/ayuth/install-anaconda-on-macos-with-homebrew-c94437d63a37

	c) fix permission on .conda file: (substitute $USER for your username)
		sudo chown -R $USER ~/.conda
	
		
2) Setup environment.

	a) cd to "TheZeroBetaRate" directory

	b) create a virtual environment with bio-env.md :
	
		conda create --name zbenv --file requirements.md
	
	c) activate virtual environment
	
		conda activate zbenv

	d) Install additional (non-conda) packages

		pip install wrds
		pip install pyarrow
	
3) Run code

	a) cd to DataCode
	
	b) edit path_variables.md to set:
		i) main_path to the "TheZeroBetaRate" directory. For example,
			 main_path: "/Users/bhebert/Documents/GitHub/TheZeroBetaRate"
		ii) your WRDS username (in quotes)
	
	c) python3 A1_CCM_download.py 
		i) on this first run, you will be asked to create a ".pgpass" file. the program will work regardless of whether you do this, but doing so will eliminate the need to enter a password each time.
		ii) on the first run, or if you didn't create a .pgpass, you will need to enter your WRDS password.
		iii) running this takes from several up to 20m in our experience, depending on your internet bandwidth to the WRDS servers
	
	d) python3 A2_beta_estimation.py
	
		i) you may get the warning "UserWarning: A worker stopped while some jobs were given to the executor." which can be disregarded. This occurs because the code tries to use all available cores to speed up this process.
	
	e) python3 00_create_data.py

		i) you may get the warning "Workbook contains no default style, apply openpyxl's default"" which can be disregarded.