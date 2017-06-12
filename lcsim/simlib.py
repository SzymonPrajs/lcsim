"""
Load, parse and store SIMLIB data in a cache file stored by default in
$HOME/.lcsim/lcsim.db
"""
import os
import sqlite3
import hashlib
import pandas as pd
import numpy as np
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

    Parameters
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

    reload_simlib : bool, optional
        Flag specifying if the SIMLIB file should be fully reloaded into the
        database even if it is already been found in the lcfit.db file.
    """
    def __init__(self, simlib_file=None, reload_simlib=False):
        home_dir = os.path.expanduser('~')
        if not os.path.exists(home_dir + '/.lcsim'):
            os.mkdir(home_dir + '/.lcsim')

        lcsim_db_path = home_dir + '/.lcsim/lcsim.db'
        self.conn = sqlite3.connect(lcsim_db_path)
        self.cur = self.conn.cursor()

        self.sqlalchemy_engine = sq.create_engine('sqlite:////'+lcsim_db_path)

        self.survey = "None"
        self._ccd_field = [0, 'C0']
        self._ra = -999
        self._dec = -999
        self._mwebv = 0.0
        self._source = [0, 0, 'n', 0, 0, 0, 0, 0, 0, 0, 0]
        self._template = {}

        self.verify_checksum_table()
        if simlib_file is not None:
            self.load_simlib_file(simlib_file, reload_simlib)

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
            table_name varchar(10),
            survey varchar(10)
        )"""
        self.cur.execute(query)
        self.conn.commit()

    def initialise_temp_simlib_table(self):
        """
        Setup a SIMLIB table in the lcsim.db SQLite database

        Parameters
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
            skysigt float,
            ra float,
            dec float,
            mwebv float
        )"""
        self.cur.execute(query)
        self.conn.commit()

    def load_simlib_file(self, simlib_file, reload_simlib=False):
        """
        Read and load a SIMLIB file into the lcsim.db database file. File is
        first MD5 hashed and looked up in the database file. If exists it is
        loaded unless reload=True in which case it is reloaded fully as in the
        case where the file has not been found before.


        Parameters
        ----------
        simlib_file : string
            Path to the SIMLIB file to be loaded.

        reload : bool, optional
            Flag specifying if the table should be fully reloaded even if it is
            already found in the lcfit.db file.

        Returns
        -------
        None
        """
        checksum = hashlib.md5(open(simlib_file, 'rb').read()).hexdigest()

        if reload_simlib is True:
            query = """
            DELETE FROM checksum where hash='{}'
            """
            self.cur.execute(query.format(checksum))
            self.conn.commit()

        query = """\
        SELECT * FROM checksum WHERE hash='{}'
        """
        self.cur.execute(query.format(checksum))
        res = self.cur.fetchall()

        if len(res) != 0:
            self.simlib_table = res[0][2]
            self.survey = res[0][3]

        else:
            self.simlib_table = 'simlib_' + checksum[0:10]
            self.initialise_temp_simlib_table()

            query = """\
            INSERT INTO checksum values{}
            """
            insert_values = tuple([checksum,
                                   simlib_file,
                                   self.simlib_table,
                                   self.survey])
            self.cur.execute(query.format(insert_values))
            self.conn.commit()

            self.parse_simlib_file(simlib_file)
            self.create_current_simlib_table()

    def update_input_values(self, key, value):
        """
        Based on the input key and value update the class wide variables that
        will eventually be pushed to the lcsim.db database file

        Parameters
        ----------
        key : str
            SNANA SIMLIB key value
        value : array-like
            Array of values corresponding to a given key
        """
        if key == 'SURVEY':
            self.survey = value[0]
            query = """\
            UPDATE checksum SET survey='{}' WHERE table_name='{}'
            """
            self.cur.execute(query.format(self.survey, self.simlib_table))
            self.conn.commit()

            for flt in value[2]:
                self._template[flt] = [0, 0]

        elif key == 'FIELD':
            self._ccd_field[1] = value[0]

            if len(value) > 3:
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

        elif key == 'RA':
            self._ra = value[0]
            self._dec = value[2]

            if len(value) > 6:
                self._mwebv = value[6]

        elif key == 'S':
            self._source = value

            if self.survey == 'DES':
                self._source[6] = np.round(float(self._source[6])*np.pi/2, 2)

    def assign_ccd_from_ra_dec(self):
        """
        Assign a CCD value for a unique pair of RA and DEC in the SIMLIB file.
        This is essentilly used as an index for the obslog for SDSS.
        """
        for field in self.get_fields():
            ccd = 1
            for pair in self.get_ccd_ra_dec(field).values:
                query = """
                UPDATE {} SET ccd={} WHERE field='{}' AND ra={} AND dec={}
                """
                self.cur.execute(query.format(self.simlib_table,
                                              ccd,
                                              field,
                                              pair[1],
                                              pair[2]))
                ccd += 1

            self.conn.commit()

    def create_current_simlib_table(self):
        """
        Depending on which survey is being used either remove duplicated
        entried (DES) or rename the 'temp' table to the correct simlim_CHECKSUM
        name in the lcsim.db. Tables are also indexed here to improve query
        performance. Finally temporary tables are removed.
        """

        query = """
        DROP TABLE IF EXISTS {}
        """
        self.cur.execute(query.format(self.simlib_table))
        self.conn.commit()

        if self.survey == 'DES':
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
                    max(skysigt) AS skysigt,
                    sum(ra) / count(ra) AS ra,
                    sum(dec) / count(dec) AS dec,
                    max(mwebv) as mwebv
                FROM temp
                GROUP BY ccd, field, idexpt, flt
            """

        else:
            query = """
            ALTER TABLE temp RENAME TO {}
            """

        self.cur.execute(query.format(self.simlib_table))
        self.conn.commit()

        index = """
        CREATE INDEX ra_dec_{}
            ON {}(ra, dec)
        """
        self.cur.execute(index.format(self.simlib_table, self.simlib_table))
        self.conn.commit()

        if self.survey == 'SDSS':
            self.assign_ccd_from_ra_dec()

        index = """
        CREATE INDEX search_query_{}
            ON {}(field, ccd, mjd, flt)
        """
        self.cur.execute(index.format(self.simlib_table, self.simlib_table))
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
            zps,sigzps,mag,zpt,skysigt,ra,dec,mwebv
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
                        tuple(self._template[self._source[2]]) +
                        (self._ra, self._dec, self._mwebv)
                    )
                    self.cur.execute(query.format(insert_values))

        self.conn.commit()

    def get_fields(self):
        """
        Get a list of unique fields available in the SIMLIB file

        Parameters
        ----------
        None

        Returns
        -------
        fields : numpy.array
            Array of unique fields found in the SIMLIB file
        """
        query = """\
            SELECT DISTINCT field
            FROM {}
            ORDER BY field asc
        """
        query = query.format(self.simlib_table)

        res = pd.read_sql_query(query, self.sqlalchemy_engine)
        res = res[res['field'] != '']
        return res['field'].values

    def get_ccds(self, field):
        """
        Get a list of unique CCDs available for a given field in the
        SIMLIB file.

        Parameters
        ----------
        field : str
            Observing field to be searched for available CCDs

        Returns
        -------
        ccds : numpy.array
            Array of unique CCDs found for the input field in the SIMLIB file
        """
        query = """\
            SELECT DISTINCT ccd
            FROM {}
            WHERE field = '{}'
            ORDER BY ccd
        """
        query = query.format(self.simlib_table, field)

        res = pd.read_sql_query(query, self.sqlalchemy_engine)
        if not np.issubdtype(res['ccd'].dtype, np.number):
            res = res[res['ccd'] != '']

        return res['ccd'].values

    def get_ccd_ra_dec(self, field):
        """
        Get the RA and DEC associated with each CCD number. For DES this is
        approximately the central RA and DEC of each CCD while for SDSS, where
        the CCD number are made up values used as indexes for RA and DEC, this
        is a usuful look up table for the number and range of available CCD
        values.

        Parameters
        ----------
        field : str
            observing field to be searched for available CCDs

        Returns
        -------
        ccd_ra_dec : Pandas.DataFrame
            DataFrame object containing unique CCD numbers and their
            corresponding RA and DEC values.
        """
        query = """\
            SELECT DISTINCT ccd, ra, dec
            FROM {}
            WHERE field = '{}'
            ORDER BY ccd, ra asc
        """
        query = query.format(self.simlib_table, field)
        return pd.read_sql_query(query, self.sqlalchemy_engine)

    def get_obslog(self, field, ccd, band=None, min_mjd=None, max_mjd=None):
        """
        Get the observing log from the cached SIMLIB file

        Parameters
        ----------
        field : str
            Name of the observing field

        ccd : int
            CCD number, must be provided as each CCD has a different set of
            observation log parameters.

        band : str or array-like, optional
            Names of observing filters. If specified the observing log will be
            returned only for these bands.

        min_mjd : float, optional
            Lower limit for the returned MJD, can be used to create light
            curves that are shorter than the full duration of the survey

        max_mjd : float, optional
            Max limit for the returned MJD, can be used to create light
            curves that are shorter than the full duration of the survey

        Returns
        -------
        obslog : Pandas.DataFrame
            `Pandas.DataFrame` returning the observation logs. The returned
            object contains the survey nama as a metadata object that can be
            accessed though obslog.survey
        """
        query = """\
            SELECT * FROM {} WHERE field="{}" AND ccd={} {} {} {}
        """
        band_str = ""
        min_mjd_str = ""
        max_mjd_str = ""

        if band is not None:
            if type(band) == str:
                band = [band]

            band_str = 'AND flt in ('
            for flt in band:
                band_str = band_str + '"' + str(flt) + '",'
            band_str = band_str[:-1] + ')'

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

        obslog = pd.read_sql_query(query, self.sqlalchemy_engine)
        obslog.survey = self.survey

        return obslog
