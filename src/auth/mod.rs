//! Secure authentication and credential management

pub mod credentials;

pub use credentials::{aws_profile, CredentialSource, Credentials};
