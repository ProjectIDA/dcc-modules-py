from ida.db import DB_TYPES

def read(db_type, db_table):

    table_data = None

    db_type = db_type.lower()

    if db_type not in DB_TYPES:
        raise ValueError('Invalid DB_TYPE: ' + db_type)

    if db_type == 'datascope': 
        import ida.db.datascope.io
        table_data = ida.db.datascope.io.read(db_table)

    return table_data



