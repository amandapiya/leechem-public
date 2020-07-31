def count_rows(table, db, before=True):
    # helper function printing number of rows in a table
    if before: text='Before'
    else: text='After'
    c = db.db.cursor()
    c.execute(f'select * from {table}')
    print(table, text, len(c.fetchall()))
