"""
Load, parse and store SIMLIB data in a cache file stored by default in
$HOME/.lcsim/lcsim.db
"""
import os
import sqlite3
import hashlib
import pandas as pd
import sqlalchemy as sq


def split_textfile_line(text_line):
    """
    Clear a line of text from spaces and linebreak characters the split to list
    """
    out = text_line.rstrip('\n')
    while '  ' in out:
        out = out.replace('  ', ' ')
    out = out.split(' ')

    return out[0][:-1], out[1:]


class SIMLIBReader():
    """
    Read SNANA's SIMLIB files and store their contents in a SQLite database

    Properties
    ----------
    simlib_file : str, optional
        Path to a SNANA SIMLIB file. If the file has already been cached
        loading will be a very quick. An md5 checksum is calculated for the
        input SIMLIB file and compared with that stored in the ~/.lcsim
        directory. If the file does not match a previously cashed version
        it will be parsed. This takes several minutes and could be heavily
        optimised using Numpa or Cython as a "one of" task the added level of
        complexity would be an unneccessary sacrifice.

        If the SIMLIB file is not passed to the constructor it can be parsed at
        any point using the `load_simlib_file(simlib_file)` method.
    """
    def __init__(self, simlib_file=None):
        home_dir = os.path.expanduser('~')
        if not os.path.exists(home_dir + '/.lcsim'):
            os.mkdir(home_dir + '/.lcsim')

        lcsim_db_path = home_dir + '/.lcsim/lcsim.db'
        self.conn = sqlite3.connect(lcsim_db_path)
        self.cur = self.conn.cursor()

        self.sqlalchemy_engine = sq.create_engine('sqlite:////'+lcsim_db_path)

        self._ccd_field = [0, 'C0']
        self._source = [0, 0, 'n', 0, 0, 0, 0, 0, 0, 0, 0]
        self._template = {'g': [0, 0], 'r': [0, 0], 'i': [0, 0], 'z': [0, 0]}

        self.verify_checksum_table()
        if simlib_file is not None:
            self.load_simlib_file(simlib_file)

    def __del__(self):
        self.conn.close()

    def verify_checksum_table(self):
        """
        Verify or create a table of checksums for SIMLIB files in the lcsim.db
        """
        query = """\
        CREATE TABLE IF NOT EXISTS checksum(
            hash varchar(32),
            file varchar(200),
            table_name varchar(10)
        )"""
        self.cur.execute(query)
        self.conn.commit()

    def initialise_temp_simlib_table(self):
        """
        Setup a SIMLIB table in the lcsim.db SQLite database

        Properties
        ----------
        table_name : str
            Name of the new SIMLIB data table

        Returns
        -------
        None
        """
        self.cur.execute("DROP TABLE IF EXISTS temp")
        self.conn.commit()

        query = """\
        CREATE TABLE IF NOT EXISTS temp (
            ccd int,
            field varchar(2),
            mjd float,
            idexpt int,
            flt varchar(1),
            gain float,
            noise float,
            skysigs float,
            psf1 float,
            psf2 float,
            psfratio float,
            zps float,
            sigzps float,
            mag float,
            zpt float,
            skysigt float
        )"""
        self.cur.execute(query)
        self.conn.commit()

    def load_simlib_file(self, simlib_file):
        """
        Read and laod a SIMLIB file into the lcsim.db database file
        """
        checksum = hashlib.md5(open(simlib_file, 'rb').read()).hexdigest()

        query = """\
        SELECT * FROM checksum WHERE hash='{}'
        """
        self.cur.execute(query.format(checksum))
        res = self.cur.fetchall()

        if len(res) != 0:
            self.simlib_table = res[0][2]

        else:
            self.simlib_table = 'simlib_' + checksum[0:10]
            self.initialise_temp_simlib_table()

            query = """\
            INSERT INTO checksum values{}
            """
            insert_values = tuple([checksum, simlib_file, self.simlib_table])
            self.cur.execute(query.format(insert_values))
            self.conn.commit()

            self.parse_simlib_file(simlib_file)
            self.create_current_simlib_table()

    def update_input_values(self, key, value):
        """
        Based on the input key and value update the class wide variables that
        will eventually be pushed to the lcsim.db database file

        Properties
        ----------
        key : str
            SNANA SIMLIB key value
        value : array-like
            Array of values corresponding to a given key
        """
        if key == 'FIELD':
            self._ccd_field[1] = value[0]
            self._ccd_field[0] = value[3][1:-1]

        elif key == 'TEMPLATE_ZPT':
            self._template['g'][0] = value[0]
            self._template['r'][0] = value[1]
            self._template['i'][0] = value[2]
            self._template['z'][0] = value[3]

        elif key == 'TEMPLATE_SKYSIG':
            self._template['g'][1] = value[0]
            self._template['r'][1] = value[1]
            self._template['i'][1] = value[2]
            self._template['z'][1] = value[3]

        elif key == 'S':
            self._source = value

    def create_current_simlib_table(self):
        """
        Remove duplicated entried in thhe lcsim.db table
        """
        query = """\
        CREATE TABLE IF NOT EXISTS {} AS
            SELECT
                ccd,
                field,
                max(mjd) AS mjd,
                idexpt,
                flt,
                max(gain) AS gain,
                max(noise) AS noise,
                max(skysigs) AS skysigs,
                max(psf1) AS psf1,
                max(psf2) AS psf2,
                max(psfratio) AS psfratio,
                max(zps) AS zps,
                max(sigzps) AS sigzps,
                max(zpt) AS zpt,
                max(skysigt) AS skysigt
            FROM temp
            GROUP BY ccd, field, idexpt, flt
        """
        self.cur.execute(query.format(self.simlib_table))
        self.conn.commit()

        self.cur.execute("DROP TABLE IF EXISTS temp")
        self.conn.commit()

    def parse_simlib_file(self, simlib_file):
        """
        Parse each line of a SIMLIB files and push to the database
        """

        query = """\
        INSERT INTO temp (
            ccd,field,mjd,idexpt,flt,gain,noise,skysigs,psf1,psf2,psfratio,
            zps,sigzps,mag,zpt,skysigt
        ) VALUES{}
        """

        with open(simlib_file, 'r') as simlib:
            for simlib_line in simlib:
                if simlib_line[0] == '#':
                    continue

                simlib_key, simlib_values = split_textfile_line(simlib_line)
                self.update_input_values(simlib_key, simlib_values)

                if simlib_key == 'S':
                    insert_values = (
                        tuple(self._ccd_field) +
                        tuple(self._source) +
                        tuple(self._template[self._source[2]])
                    )
                    self.cur.execute(query.format(insert_values))

        self.conn.commit()

    def get_obslog(self, field, ccd, band=None, min_mjd=None, max_mjd=None):
        """
        Get the observing log from the cached SIMLIB file

        Properties
        ----------
        field : str
            Name of the observing field

        ccd : int
            CCD number, must be provided as each CCD has a different set of
            observation log parameters.

        band : str, optional
            Name of an observed filter. If specified the observing log will be
            returned only for this band.

        min_mjd : float, optional
            Lower limit for the returned MJD, can be used to create light
            curves that are shorter than the full duration of the survey

        max_mjd : float, optional
            Max limit for the returned MJD, can be used to create light
            curves that are shorter than the full duration of the survey

        Returns
        -------
        obslog : Pandas.DataFrame
            `Pandas.DataFrame` returning the observation logs.
        """
        query = """\
            SELECT * FROM {} WHERE field="{}" AND ccd={} {} {} {}
        """
        band_str = ""
        min_mjd_str = ""
        max_mjd_str = ""

        if band is not None:
            band_str = 'AND flt="' + str(band) + '"'

        if min_mjd is not None:
            min_mjd_str = 'AND mjd>' + str(min_mjd)

        if max_mjd is not None:
            max_mjd_str = 'AND mjd<' + str(max_mjd)

        query = query.format(self.simlib_table,
                             field,
                             ccd,
                             band_str,
                             min_mjd_str,
                             max_mjd_str)
        return pd.read_sql_query(query, self.sqlalchemy_engine)
