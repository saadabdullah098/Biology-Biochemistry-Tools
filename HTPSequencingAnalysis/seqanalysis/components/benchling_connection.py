import os

from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth

class BenchlingConnection:
    """
    Parameters
    ----------
    benchling api password and url that can be called as variables

    Results
    -------
    returns api connection string
    """

    def __init__(self, environment: str = None):
        self.environment = environment

    def benchling_connection(self):
        """
        calls connection variables based on user input
        """
        if "production" == self.environment:
            benchling = self.prod_con()

        elif "test" == self.environment:
            benchling = self.test_con()

        else:
            raise Exception("environment input must be 'production' or 'test' and is case sensitive")

        return benchling

    def test_con(self):
        """
        sets test connection string
        """
        url = os.getenv("test_url")
        api_auth = os.getenv("test_api_key")
        benchling = Benchling(url=url, auth_method=ApiKeyAuth(api_auth))

        return benchling

    def prod_con(self):
        """
        sets production connection string
        """
        url = os.getenv("prod_url")
        api_auth = os.getenv("prod_api_key")
        benchling = Benchling(url=url, auth_method=ApiKeyAuth(api_auth))

        return benchling
