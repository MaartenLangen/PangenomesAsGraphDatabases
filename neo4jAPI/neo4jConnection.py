import neo4j

class Connection:
    
    def __init__(self, uri: str, user: str, pwd: str):
        """Initialization of a class' instance. Sets uri, user, and pwd and
        creates a connection to the database.
        
        ## Input
        - uri: uri to the database
        - user: username for authentication
        - pwd: password for authentication"""

        self.__uri = uri
        self.__user = user
        self.__pwd = pwd
        self.__driver = None

        try:
            self.__driver = neo4j.GraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
        except Exception as e:
            print("Failed to create the driver:", e)
        
    def close(self):
        """Close the driver."""

        if self.__driver is not None:
            self.__driver.close()

    def getDriver(self) -> neo4j.Driver:
        """Retrieve the driver.
        
        ## Output: 
        - the driver of this instance"""

        return self.__driver