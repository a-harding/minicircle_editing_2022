"""Functions for sending email notifications with updates about program stage completion, preliminary
analysis, premature termination errors, and completion notification."""


import os
import smtplib
import ssl
from email.message import EmailMessage
from email.headerregistry import Address


def email_notification(message: str, type: str=''):
    """Sends an email containing the supplied message string."""

    port = 465 # for use with smtplib.SMTP_SSL()
    recipient = 'adam0harding@gmail.com'

    message1 = """
    Subject: Subject text.
    
    Some text."""

    email_deet_path = '/home/adam/Documents/SchnauferLab/Rotation1/minicircle_editing/minicircle_editing/support_files/coding_email.txt'
    with open(email_deet_path, 'r') as f:
        eaddress, epass = f.read().split()

    msg = EmailMessage()
    msg['Subject'] = 'Subject text.'
    msg['From'] = eaddress
    msg['To'] = recipient
    msg.set_content('Content text')

    #  CREATES A SSL-SECURED CONTEXT
    context = ssl.create_default_context()

    with smtplib.SMTP_SSL('smtp.gmail.com', port=port, context=context) as server:
        server.login(eaddress, epass)
        server.sendmail(from_addr=eaddress, to_addrs=recipient, msg=message1)