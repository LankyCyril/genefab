def safe_file_name(url):
    return url.split("/")[-1].lstrip("~")
