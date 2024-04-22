import requests


def get_info():
    res = requests.get('http://127.0.0.1:8000/api/')
    if res.status_code == 200:
        return res.json()
    else:
        return None
