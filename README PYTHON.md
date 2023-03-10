# TheZeroBetaRate


Python-related instructions for Mac OS X:

1) Install python. One set of reasonable instructions:

	a) install homebrew
	
	b) install anaconda with homebrew
		https://medium.com/ayuth/install-anaconda-on-macos-with-homebrew-c94437d63a37

	c) fix permission on .conda file: (substitute $USER for your username)
		sudo chown -R $USER ~/.conda
	
	d) create a virtual environment:
	
		conda create --name zbenv
		
	e) activate virtual environment
	
		conda activate zbenv
		
		
2) Install packages

	conda install pandas
	conda install xlrd
	conda install openpyxl
	pip install wrds
	pip install pyarrow


3) Run code

	a) cd to "TheZeroBetaRate" directory
	
	b) cd to DataCode
	
	c) edit 00_create_data.py to set:
		i) main_path to the "TheZeroBetaRate" directory
		ii) your WRDS username
	
	d) python3 00_create_data.py
	
		i) on this first run, you will be asked to create a ".pgpass" file. the program will work regardless of whether you do this, but doing so will eliminate the need to enter a password each time.
		ii) on the first run, or if you didn't create a .pgpass, you will need to enter your WRDS password.
		iii) you will get the warning "Workbook contains no default style, apply openpyxl's default"" which can be disregarded.
	