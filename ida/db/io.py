from ida.db import DB_TYPES

def read(db_type, db_table):

    table_data = None
    result = False

    db_type = db_type.lower()

    if db_type not in DB_TYPES:
        raise ValueError('Invalid DB_TYPE: ' + db_type)
    elif db_type == 'datascope':
        import ida.db.datascope.io
        result, table_data = ida.db.datascope.io.read(db_table)

    return result, table_data



