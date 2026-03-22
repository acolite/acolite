//! Secure authentication and credential management

pub mod credentials;

pub use credentials::{aws_profile, CdseCredentials, CredentialSource, Credentials};
