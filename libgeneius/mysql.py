from libgeneius.error import GeneiusError
try:
    import MySQLdb
except:
    raise GeneiusError("Missing Package: MySQLdb")

class GeneiusDb:

    def __init__(self,mydb,myhost,myuser,mypasswd):
                
        try:
            self.db = MySQLdb.connect(host=myhost,
                                 user=myuser,
                                 passwd=mypasswd,
                                 db = mydb)
        except MySQLdb.OperationalError, oe:
            raise GeneiusError("Database error %s" % str(oe))
        
        self.cursor = self.db.cursor()

    def query(self,qstring):
        self.cursor.execute(qstring)
        return self.cursor.fetchall()
