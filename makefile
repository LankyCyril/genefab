genefab/_readme.py: assets/readme.html
	echo 'html = """' > $@
	cat $< >> $@
	echo '"""' >> $@
