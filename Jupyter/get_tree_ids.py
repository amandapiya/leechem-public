def get_tree_ids(db, where=''):
    '''Gets a tuple of tree IDs that are in db.'''

    query = "select tid from trees"

    if where != '':
        query += f' where {where}'

    c = db.db.cursor()
    c.execute(query)
    all_rows = c.fetchall()
    return list(list(zip(*all_rows))[0])
