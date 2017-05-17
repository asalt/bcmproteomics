from concurrent.futures import ThreadPoolExecutor

from . import ispec

def get_data(recno=None, runno=1, searchno=1, typeof=ispec.E2G, verbosity=1, **kwargs):
    """asyncronous data access for an experiment"""
    if not any(typeof==x for x in (ispec.Experiment, ispec.E2G, ispec.PSMs)):
        print('Invalid input `typeof`', typeof)
        return None
    if verbosity > 0:
        print('Getting {} data for {}_{}_{}'.format(typeof.__name__,
                                                    recno,
                                                    runno,
                                                    searchno))
    return typeof(recno, runno, searchno, **kwargs)

def fmt_input_data(datas,):
    """datas is an iterable of rec, run, search
    """
    if isinstance(datas, (str, int, float)):
        raise ValueError('input must be list-like')
    fmt_datas = list()
    for data in datas:
        if isinstance(data, dict):
            fmt_datas.append(data)
        elif isinstance(data, tuple):
            try:
                recno, runno, searchno = data
            except ValueError:
                print('Invalid format for', data)
                continue
            fmt_datas.append(dict(recno=recno,
                                  runno=runno,
                                  searchno=searchno))
        elif isinstance(data, (int, float, str)):
            fmt_datas.append(dict(recno=data))
    return fmt_datas

def async_get_datas(datas, verbosity=1, typeof=None, max_workers=10, **kwargs):
    """
    :datas: an iterable of rec, run, search
    typeof: one of ispec.Experiment, ispec.E2G, or ispec.PSMs
    """
    if typeof is None:
        typeof = ispec.E2G
    if not any(typeof==x for x in (ispec.Experiment, ispec.E2G, ispec.PSMs)):
        raise TypeError('Invalid input `typeof`', typeof)

    conn = ispec.filedb_connect()
    if isinstance(conn, str):
        raise Exception(conn)

    fmt_data = fmt_input_data(datas)
    if not fmt_data:
        return

    futures = list()
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(get_data, **data, verbosity=verbosity, typeof=typeof, **kwargs)
                   for data in fmt_data]

    return [future.result() for future in futures]

    # loop = asyncio.get_event_loop()
    # loop.run_until_complete

    # coros = [run_get_data(**data, verbosity=verbosity, typeof=typeof, **kwargs)
    #          for data in fmt_data]
    # fg = asyncio.gather( *coros )
    # results = loop.run_until_complete(fg)
    # return results
